#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Somatic V/D/J recombination analyzer
# Author: Lior Galanti < lior.galanti@nyu.edu >
# NYU Center for Genetics and System Biology 2015
#     
# somatic is free software; you can redistribute it and/or modify it under the terms of 
# the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.
# 

import sys
import logging
import re
import json
import uuid
import io
import hashlib 

from io import StringIO, BytesIO
from datetime import timedelta, datetime
from argparse import ArgumentParser
from subprocess import Popen, PIPE

from bson.objectid import ObjectId
from pymongo.son_manipulator import SONManipulator
from pymongo import MongoClient, DESCENDING, ASCENDING
from pymongo.errors import BulkWriteError

def to_json(node):
    def handler(o):
        result = None
        if isinstance(o, datetime):
            result = o.isoformat()
        if isinstance(o, ObjectId):
            result = str(o)
        if isinstance(o, set):
            result = list(o)
        if isinstance(o, Sequence):
            result = o.node
        if isinstance(o, Sample):
            result = o.head
        return result

    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4, default=handler)

log_levels = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

configuration = {
    'diagram': {
        'prototype': {
            'name': {
                'format': lambda x: x,
                'title': 'name',
                'width': 'auto',
                'value': 'subject id'
            },
            'strain': {
                'format': lambda x: x,
                'title': 'strain',
                'width': 'auto',
                'value': 'strain',
            },
            'region': {
                'format': lambda x: x,
                'title': 'R',
                'width': 'auto',
                'value': 'region',
            },
            'functionality': {
                'format': lambda x: x,
                'title': 'F',
                'width': 1,
                'value': 'functionality',
            },
            'in frame': {
                'format': lambda x: '*' if x else ' ',
                'title': 'I',
                'width': 1,
                'value': 'in frame',
            },
            'picked': {
                'format': lambda x: '*' if x else ' ',
                'title': 'P',
                'width': 1,
                'value': 'picked',
            },
            'gapped': {
                'format': lambda x: '*' if x else ' ',
                'title': 'G',
                'width': 1,
                'value': 'gapped',
            },
            'strand': {
                'format': lambda x: '+' if x else '-',
                'title': 'S',
                'width': 1,
                'value': 'subject strand',
            }            
        }
    },
    'profile': {
        'b6~': {
            'reference': {
                'head.verified': True,
                'head.organism name': 'mus musculus',
                '$or': [ { 'head.strain': { '$ne': 'C57BL/6' } }, { 'head.strain': { '$exists': False } } ]
            },
        },
        'b6': {
            'reference': {
                'head.verified': True,
                'head.organism name': 'mus musculus',
                'head.strain': 'C57BL/6'
            },
        },
        'mouse': {
            'reference': {
                'head.verified': True,
                'head.organism name': 'mus musculus',
            },
        },
        'default': {
            'reference': {
                'head.verified': True,
            },
            'sample': {
                'head.valid': True,
                'head.productive': True,
            },
            'diagram': {
                'track': {
                    'valid': True
                },
                'feature': [
                    'region',
                    'name',
                    'strain',
                    'strand',
                    'functionality',
                    'in frame',
                    'gapped',
                    'picked',
                ]
            },
        },
        'productive': {
            'sample': {
                'head.valid': True,
                'head.productive': True,
            },
            'diagram': {
                'track': {
                    'valid': True,
                    'picked': True
                },
                'feature': [
                    'region',
                    'name',
                    'strain',
                    'functionality',
                    'strand',
                ]
            },
        },
        'productive!': {
            'sample': {
                'head.valid': True,
                'head.productive': False,
            },
            'diagram': {
                'track': {
                    'valid': True,
                    'picked': True
                },
                'feature': [
                    'region',
                    'name',
                    'strain',
                    'strand',
                    'functionality',
                    'in frame',
                ]
            },
        },
        'productive~': {
            'sample': {
                'head.valid': True,
                'head.productive': True,
            },
            'diagram': {
                'track': {
                    'valid': True,
                    'picked': False
                },
                'feature': [
                    'region',
                    'name',
                    'strain',
                    'strand',
                    'functionality',
                    'in frame',
                ]
            },
        },
        'palindromic': {
            'sample': {
                'head.valid': True,
                'head.productive': True,
                'head.palindromic': True,
            },
            'diagram': {
                'track': {
                    'valid': True,
                    'picked': True
                },
                'feature': [
                    'region',
                    'name',
                    'strain',
                    'functionality',
                    'strand',
                ]
            },
        },
        'gapped': {
            'sample': {
                'head.valid': True,
                'head.gapped': True,
            },
            'diagram': {
                'track': {
                    'valid': False
                },
                'feature': [
                    'region',
                    'name',
                    'strain',
                    'functionality',
                    'strand',
                ]
            },
        },
        'premature': {
            'sample': {
                'head.valid': True,
                'head.premature': True,
            },
            'diagram': {
                'track': {
                    'valid': True,
                    'picked': True
                },
                'feature': [
                    'region',
                    'name',
                    'strain',
                    'strand',
                    'functionality',
                    'in frame',
                ]
            },
        },
        'in': {
            'sample': {
                'head.valid': True,
                'head.in frame': True,
            },
            'diagram': {
                'track': {
                    'valid': True,
                    'picked': True
                },
                'feature': [
                    'region',
                    'name',
                    'strain',
                    'strand',
                    'functionality',
                ]
            },
        },
        'out': {
            'sample': {
                'head.valid': True,
                'head.in frame': False,
            },
            'diagram': {
                'track': {
                    'valid': True,
                    'picked': True
                },
                'feature': [
                    'region',
                    'name',
                    'strain',
                    'strand',
                    'functionality',
                ]
            },
        },
    },
    'yes/no question': ['Y', 'N'],
    'regions': ['VH', 'DH', 'JH'],
    'expression': {
        'gapped sequence': re.compile('^\s+Query_[0-9]+\s+(?P<offset>[0-9]+)\s+(?P<sequence>[ATCGN-]+)\s+[0-9]+$'),
        'expand hit': re.compile(
            r"""
            (?P<region>VH|DH|JH),
            (?P<subject_id>[^,]+),
            (?P<query_start>[^,]+),
            (?P<query_end>[^,]+),
            (?P<subject_start>[^,]+),
            (?P<subject_end>[^,]+),
            (?P<gap_openings>[^,]+),
            (?P<gaps>[^,]+),
            (?P<mismatch>[^,]+),
            (?P<identical>[^,]+),
            (?P<bit_score>[^,]+),
            (?P<evalue>[^,]+),
            (?P<alignment_length>[^,]+),
            (?P<subject_strand>[^,]),
            (?P<query_strand>[^,])
            """,
            re.VERBOSE
        ),
        'igblast compressed hit': '{region},{subject id},{query start},{query end},{subject start},{subject end},{gap openings},{gaps},{mismatch},{identical},{bit score},{evalue},{alignment length},{subject strand},{query strand}',
        'igblast hit': re.compile(
            r"""
            (?P<region>[VDJ])\t
            (?:reversed\|)?(?:[^,]+)\t
            (?P<subject_id>[^,]+)\t
            (?P<query_start>[^,]+)\t
            (?P<query_end>[^,]+)\t
            (?P<subject_start>[^,]+)\t
            (?P<subject_end>[^,]+)\t
            (?P<gap_openings>[^,]+)\t
            (?P<gaps>[^,]+)\t
            (?P<mismatch>[^,]+)\t
            (?P<identical>[^,]+)\t
            (?P<bit_score>[^,]+)\t
            (?P<evalue>[^,]+)\t
            (?P<alignment_length>[^,]+)\t
            (?P<subject_strand>[^,]+)
            """,
            re.VERBOSE
        ),
        'igblast reversed query': '# Note that your query represents the minus strand of a V gene and has been converted to the plus strand.',
        'imgt fasta header': re.compile(
            r"""
            ^>
            (?P<accession>[^|]*)\|
            (?P<allele_name>(?P<gene_name>(?P<subgroup>[^-]+)[^\*]+)\*[0-9]+)\|
            (?P<organism_name>[^|_]*)(?:_(?P<strain>[^|]+))?\|
            (?P<functionality>[\[\(]?(?:F|ORF|P)[\]\)]?)\|
            (?P<region>V|D|J)-REGION\|
            (?:(?P<start>[0-9]+)\.\.(?P<end>[0-9]+))?\s*\|
            (?P<length>[0-9]+)\snt\|
            (?P<read_frame>[123]|NR)\|
            [^|]*\|
            [^|]*\|
            [^|]*\|
            [^|]*\|
            [^|]*\|
            [^|]*\|
            (?P<polarity>rev-compl)?\s*\|$
            """,
            re.VERBOSE
        ),
        'nucleotide sequence': re.compile('^[ACGTRYKMSWBDHVN]+$', re.IGNORECASE)
    },
    'accession header': {
        'imgt': re.compile(r'^>(?P<accession_number>[^ ]+) (?P<description>.*)$'),
        'ncbi': re.compile(r'^>gi\|[0-9]+\|[a-z]+\|(?P<accession_number>[^\.]+)+\.(?P<accession_version>[0-9]+)\|(?P<description>.*)$'),
    },
    'strain': {
        'mus musculus': [
            {
                'name': 'C57BL/6',
                'expression': re.compile(r'\b(?:c57bl/6|c57)\b', re.IGNORECASE)
            },
            {
                'name': 'C57BL/6J',
                'expression': re.compile(r'\bc57bl/6j\b', re.IGNORECASE)
            },
            {
                'name': 'C57BL/6N',
                'expression': re.compile(r'\bc57bl/6n\b', re.IGNORECASE)
            },
            {
                'name': 'C57BL/10',
                'expression': re.compile(r'\bc57bl/10\b', re.IGNORECASE)
            },
            {
                'name': '129/Sv',
                'expression': re.compile(r'\b129/sv\b', re.IGNORECASE)
            },
            {
                'name': '129/J',
                'expression': re.compile(r'\b129/j\b', re.IGNORECASE)
            },
            {
                'name': '129/Ole',
                'expression': re.compile(r'\b129/ole\b', re.IGNORECASE)
            },
            {
                'name': 'BALB/c',
                'expression': re.compile(r'\bbalb/c\b', re.IGNORECASE)
            },
            {
                'name': 'BALB.K',
                'expression': re.compile(r'\bbalb\.k\b', re.IGNORECASE)
            },
            {
                'name': 'MRL/lpr',
                'expression': re.compile(r'\bmrl/lpr\b', re.IGNORECASE)
            },
            {
                'name': 'A/J',
                'expression': re.compile(r'\ba/j\b', re.IGNORECASE)
            },
            {
                'name': 'NBZ',
                'expression': re.compile(r'\bnbz\b', re.IGNORECASE)
            },
        ]
    },
    'command': {
        'blat': {
            'cwd': '/Users/lg/code/blat',
            'arguments': [
                'blat',
                'chr12.fa',
                'stdin',
                'stdout',
                '-minIdentity=100'
            ]
        },
        'igblast': {
            'cwd': '/Users/lg/code/somatic/igblast',
            'arguments': [
                'igblastn',
                '-germline_db_V', 'database/mouse_imgt_vh',
                '-germline_db_J', 'database/mouse_imgt_jh',
                '-germline_db_D', 'database/mouse_imgt_dh',
                '-num_alignments_V', '7',
                '-num_alignments_J', '3',
                '-num_alignments_D', '5',
                '-organism', 'mouse',
                '-domain_system', 'imgt',
                '-query', '-',
                '-auxiliary_data', 'optional_file/mouse_gl.aux',
                '-show_translation',
                '-num_threads', '8',
                '-outfmt',
                '7 qseqid sseqid qstart qend sstart send gapopen gaps mismatch pident bitscore evalue length sstrand',
            ]
        },
        'igblast.gapped': {
            'cwd': '/Users/lg/code/somatic/igblast',
            'arguments': [
                'igblastn',
                '-germline_db_V', 'database/mouse_imgt_vh',
                '-germline_db_J', 'database/mouse_imgt_jh',
                '-germline_db_D', 'database/mouse_imgt_dh',
                '-num_alignments_V', '7',
                '-num_alignments_J', '3',
                '-num_alignments_D', '5',
                '-organism', 'mouse',
                '-domain_system', 'imgt',
                '-query', '-',
                '-auxiliary_data', 'optional_file/mouse_gl.aux',
                '-show_translation',
                '-num_threads', '8',
                '-outfmt', '4',
            ]
        } 
    },
    'interface': {
        'global': {
            'argument': [
                'version', 
                'verbosity'
            ]
        }, 
        'instruction': {
            'description': 'Lior Galanti lior.galanti@nyu.edu NYU Center for Genomics & Systems Biology'
        }, 
        'prototype': {
            'alignment': {
                'flag': [
                    '-A', 
                    '--alignment'
                ], 
                'parameter': {
                    'action': 'store_true', 
                    'dest': 'alignment', 
                    'help': 'print alignment diagram'
                }
            }, 
            'limit': {
                'flag': [
                    '-L', 
                    '--limit'
                ], 
                'parameter': {
                    'type': 'int',
                    'dest': 'limit', 
                    'help': 'max results to return'
                }
            }, 
            'skip': {
                'flag': [
                    '-S', 
                    '--skip'
                ], 
                'parameter': {
                    'type': 'int',
                    'dest': 'skip', 
                    'help': 'skip the first SKIP results'
                }
            }, 
            'gapped': {
                'flag': [
                    '-G', 
                    '--gapped'
                ], 
                'parameter': {
                    'choices': [
                        'Y', 
                        'N'
                    ], 
                    'dest': 'gapped', 
                    'help': 'gapped alignment'
                }
            }, 
            'drop': {
                'flag': [
                    '-D', 
                    '--drop'
                ], 
                'parameter': {
                    'action': 'store_true', 
                    'dest': 'drop', 
                    'help': 'drop library items'
                }
            }, 
            'strain': {
                'flag': [
                    '--strain'
                ], 
                'parameter': {
                    'dest': 'strain', 
                    'help': 'strain'
                }
            }, 
            'id': {
                'flag': [
                    '-d', 
                    '--id'
                ], 
                'parameter': {
                    'dest': 'id', 
                    'help': 'sample id', 
                    'metavar': 'ID'
                }
            }, 
            'in frame': {
                'flag': [
                    '-F', 
                    '--in-frame'
                ], 
                'parameter': {
                    'choices': [
                        'Y', 
                        'N'
                    ], 
                    'dest': 'in frame', 
                    'help': 'in frame alignment'
                }
            }, 
            'productive': {
                'flag': [
                    '-P', 
                    '--productive'
                ], 
                'parameter': {
                    'choices': [
                        'Y', 
                        'N'
                    ], 
                    'dest': 'productive', 
                    'help': 'productive alignment'
                }
            }, 
            'json': {
                'flag': [
                    '-J', 
                    '--json'
                ], 
                'parameter': {
                    'action': 'store_true', 
                    'dest': 'json', 
                    'help': 'print json info'
                }
            }, 
            'library': {
                'flag': [
                    '-l', 
                    '--library'
                ], 
                'parameter': {
                    'dest': 'library', 
                    'help': 'library name', 
                    'metavar': 'NAME'
                }
            }, 
            'format': {
                'flag': [
                    '-f', 
                    '--format'
                ], 
                'parameter': {
                    'required': True,
                    'choices': [
                        'imgt', 
                        'ncbi'
                    ], 
                    'dest': 'format', 
                    'help': 'fasta header format'
                }
            }, 
            'path': {
                'flag': [
                    'path'
                ], 
                'parameter': {
                    'help': 'file paths', 
                    'metavar': 'PATH', 
                    'nargs': '*'
                }
            }, 
            'premature': {
                'flag': [
                    '-T', 
                    '--premature'
                ], 
                'parameter': {
                    'choices': [
                        'Y', 
                        'N'
                    ], 
                    'dest': 'premature', 
                    'help': 'premature termination alignment'
                }
            }, 
            'region': {
                'flag': [
                    '-r', 
                    '--region'
                ], 
                'parameter': {
                    'choices': [
                        'VH', 
                        'DH', 
                        'JH'
                    ], 
                    'dest': 'region', 
                    'help': 'region'
                }
            }, 
            'profile': {
                'flag': [
                    '-p', 
                    '--profile'
                ], 
                'parameter': {
                    'default': 'default',
                    'choices': None,
                    'metavar': 'NAME',
                    'dest': 'profile', 
                    'help': 'profile is one of default, productive, productive.dropped, psudo, premature, inframe, outframe'
                }
            }, 
            'valid': {
                'flag': [
                    '-V', 
                    '--valid'
                ], 
                'parameter': {
                    'choices': [
                        'Y', 
                        'N'
                    ], 
                    'dest': 'valid', 
                    'help': 'valid alignment'
                }
            }, 
            'verbosity': {
                'flag': [
                    '-v', 
                    '--verbosity'
                ], 
                'parameter': {
                    'choices': [
                        'debug', 
                        'info', 
                        'warning', 
                        'error', 
                        'critical'
                    ], 
                    'default': 'info', 
                    'dest': 'verbosity', 
                    'help': 'logging verbosity level', 
                    'metavar': 'LEVEL'
                }
            }, 
            'version': {
                'flag': [
                    '--version'
                ], 
                'parameter': {
                    'action': 'version', 
                    'version': '%(prog)s 1.0'
                }
            }
        }, 
        'section': {
            'action': [
                {
                    'argument': [
                        'library',
                        'strain',
                        'drop'
                    ], 
                    'instruction': {
                        'description': 'match each read in file to regions with igblast and store results in the library. takes data from stdin', 
                        'help': 'populate samples for library', 
                        'name': 'populate'
                    }
                }, 
                {
                    'argument': [
                        'library'
                    ], 
                    'instruction': {
                        'description': 'remove all samples for the provided library name', 
                        'help': 'drop samples of a given library', 
                        'name': 'drop'
                    }
                }, 
                {
                    'argument': [
                        'profile',
                        'library', 
                        'id', 
                        'gapped', 
                        'valid', 
                        'in frame', 
                        'premature',
                        'limit',
                        'skip',
                        'productive',
                    ], 
                    'instruction': {
                        'help': 'view sample alignment', 
                        'name': 'view'
                    }
                }, 
                {
                    'argument': [
                        'profile',
                        'library', 
                        'id', 
                        'gapped', 
                        'valid', 
                        'in frame', 
                        'premature',
                        'productive',
                    ], 
                    'instruction': {
                        'help': 'count sample alignment', 
                        'name': 'count'
                    }
                }, 
                {
                    'argument': [
                        'profile',
                        'library', 
                        'id', 
                        'gapped', 
                        'valid', 
                        'in frame', 
                        'premature',
                        'limit',
                        'skip',
                        'productive',
                    ], 
                    'instruction': {
                        'help': 'print sample fasta stream', 
                        'name': 'fasta'
                    }
                }, 
                {
                    'argument': [
                        'profile',
                        'library', 
                        'id', 
                        'gapped', 
                        'valid', 
                        'in frame', 
                        'premature',
                        'limit',
                        'skip',
                        'productive',
                    ], 
                    'instruction': {
                        'help': 'print sample fastq stream', 
                        'name': 'fastq'
                    }
                }, 
                {
                    'argument': [
                        'profile',
                        'library', 
                        'id', 
                        'gapped', 
                        'valid', 
                        'in frame', 
                        'premature',
                        'limit',
                        'skip',
                        'productive',
                    ], 
                    'instruction': {
                        'help': 'view JSON sample record', 
                        'name': 'info'
                    }
                }, 
                {
                    'argument': [
                        'profile',
                        'library', 
                        'strain',
                        'json', 
                        'alignment'
                    ], 
                    'instruction': {
                        'help': 'process alignment without updating the database', 
                        'name': 'simulate'
                    }
                }, 
                {
                    'argument': [
                        'format',
                    ],
                    'instruction': {
                        'help': 'load imgt accession sequences from fasta file. takes data from stdin', 
                        'name': 'accession-populate'
                    }
                }, 
                {
                    'argument': [],
                    'instruction': {
                        'help': 'dump accession sequences to fasta', 
                        'name': 'accession-fasta'
                    }
                }, 
                {
                    'argument': [
                        'path'
                    ], 
                    'instruction': {
                        'help': 'load imgt reference sequences from fasta file', 
                        'name': 'ref-populate'
                    }
                }, 
                {
                    'argument': [
                        'region',
                        'strain',
                        'profile',
                        'id', 
                    ], 
                    'instruction': {
                        'help': 'dump reference sequences to fasta', 
                        'name': 'ref-fasta'
                    }
                }, 
                {
                    'argument': [
                        'region',
                        'strain',
                        'profile',
                        'id', 
                    ], 
                    'instruction': {
                        'help': 'dump reference sequences to igblast auxiliary file', 
                        'name': 'ref-igblast-aux'
                    }
                }, 
                {
                    'argument': [], 
                    'instruction': {
                        'help': 'rebuild database indexes', 
                        'name': 'rebuild'
                    }
                }
            ], 
            'instruction': {
                'description': '', 
                'dest': 'action', 
                'help': None, 
                'metavar': 'ACTION', 
                'title': 'pipeline operations'
            }
        }
    },
    'iupac amino acid notation': {
        'A': {  'code':'ala',   'weight': 89,   'charge': 1.8,  'name': 'Alanine' },
        'R': {  'code':'arg',   'weight': 174,  'charge': -4.5, 'name': 'Arginine' },
        'N': {  'code':'asn',   'weight': 132,  'charge': -3.5, 'name': 'Asparagine' },
        'D': {  'code':'asp',   'weight': 133,  'charge': -3.5, 'name': 'Aspartic acid' },
        'C': {  'code':'cys',   'weight': 121,  'charge': 2.5,  'name': 'Cysteine' },
        'Q': {  'code':'gln',   'weight': 146,  'charge': -3.5, 'name': 'Glutamine' },
        'E': {  'code':'glu',   'weight': 147,  'charge': -3.5, 'name': 'Glutamic acid' },
        'G': {  'code':'gly',   'weight': 75,   'charge': -0.4, 'name': 'Glycine' },
        'H': {  'code':'his',   'weight': 155,  'charge': -3.2, 'name': 'Histidine' },
        'I': {  'code':'ile',   'weight': 131,  'charge': 4.5,  'name': 'Isoleucine' },
        'L': {  'code':'leu',   'weight': 131,  'charge': 3.8,  'name': 'Leucine' },
        'K': {  'code':'lys',   'weight': 146,  'charge': -3.9, 'name': 'Lysine' },
        'M': {  'code':'met',   'weight': 149,  'charge': 1.9,  'name': 'Methionine' },
        'F': {  'code':'phe',   'weight': 165,  'charge': 2.8,  'name': 'Phenylalanine' },
        'P': {  'code':'pro',   'weight': 115,  'charge': -1.6, 'name': 'Proline' },
        'S': {  'code':'ser',   'weight': 105,  'charge': -0.8, 'name': 'Serine' },
        'T': {  'code':'thr',   'weight': 119,  'charge': -0.7, 'name': 'Threonine' },
        'W': {  'code':'trp',   'weight': 204,  'charge': -0.9, 'name': 'Tryptophan' },
        'Y': {  'code':'tyr',   'weight': 181,  'charge': -1.3, 'name': 'Tyrosine' },
        'V': {  'code':'val',   'weight': 117,  'charge': 4.2,  'name': 'Valine' },
        'B': {  'code':'asx',                                   'name': 'Aspartic acid or Asparagine' },
        'Z': {  'code':'glx',                                   'name': 'Glutamine or Glutamic acid' },
        'J': {  'code':'xle',                                   'name': 'Leucine or Isoleucine' },
        'X': {  'code':'xaa',	                                'name': 'Any amino acid' }
    },
    'iupac nucleic acid notation': {
        'A': {  'reverse': 'T',  'option':[ 'A' ],                  'name': 'Adenine' },
        'C': {  'reverse': 'G',  'option':[ 'C' ],                  'name': 'Cytosine' },
        'G': {  'reverse': 'C',  'option':[ 'G' ],                  'name': 'Guanine' },
        'T': {  'reverse': 'A',  'option':[ 'T' ],                  'name': 'Thymine' },
        'R': {  'reverse': 'Y',  'option':[ 'G', 'A' ],             'name': 'Purine' },
        'Y': {  'reverse': 'R',  'option':[ 'T', 'C' ],             'name': 'Pyrimidine' },
        'K': {  'reverse': 'M',  'option':[ 'G', 'T' ],             'name': 'Keto' },
        'M': {  'reverse': 'K',  'option':[ 'A', 'C' ],             'name': 'Amino' },
        'S': {  'reverse': 'S',  'option':[ 'G', 'C' ],             'name': 'Strong bonds' },
        'W': {  'reverse': 'W',  'option':[ 'A', 'T' ],             'name': 'Weak bonds' },
        'B': {  'reverse': 'V',  'option':[ 'G', 'T', 'C' ],        'name': 'Not A' },
        'D': {  'reverse': 'H',  'option':[ 'G', 'A', 'T' ],        'name': 'Not C' },
        'H': {  'reverse': 'D',  'option':[ 'A', 'C', 'T' ],        'name': 'Not G' },
        'V': {  'reverse': 'B',  'option':[ 'G', 'C', 'A' ],        'name': 'Not T' },
        'N': {  'reverse': 'N',  'option':[ 'A', 'G', 'C', 'T' ],   'name': 'Any' }
    },
    'nucleic to amino': {
        'GCT':'A',
        'GCC':'A',
        'GCA':'A',
        'GCG':'A',
        'TTA':'L',
        'TTG':'L',
        'CTT':'L',
        'CTC':'L',
        'CTA':'L',
        'CTG':'L',
        'CGT':'R',
        'CGC':'R',
        'CGA':'R',
        'CGG':'R',
        'AGA':'R',
        'AGG':'R',
        'AAA':'K',
        'AAG':'K',
        'AAT':'N',
        'AAC':'N',
        'ATG':'M',
        'GAT':'D',
        'GAC':'D',
        'TTT':'F',
        'TTC':'F',
        'TGT':'C',
        'TGC':'C',
        'CCT':'P',
        'CCC':'P',
        'CCA':'P',
        'CCG':'P',
        'CAA':'Q',
        'CAG':'Q',
        'TCT':'S',
        'TCC':'S',
        'TCA':'S',
        'TCG':'S',
        'AGT':'S',
        'AGC':'S',
        'GAA':'E',
        'GAG':'E',
        'ACT':'T',
        'ACC':'T',
        'ACA':'T',
        'ACG':'T',
        'GGT':'G',
        'GGC':'G',
        'GGA':'G',
        'GGG':'G',
        'TGG':'W',
        'CAT':'H',
        'CAC':'H',
        'TAT':'Y',
        'TAC':'Y',
        'ATT':'I',
        'ATC':'I',
        'ATA':'I',
        'GTT':'V',
        'GTC':'V',
        'GTA':'V',
        'GTG':'V',
        'TAA':'*',
        'TGA':'*',
        'TAG':'*',
    },
    'constant': {
        'buffer size': 32,
        'fasta line length': 80,
        'minimum v identity': 0.7, 
        'minimum d identity': 0.7, 
        'minimum j identity': 0.8, 
        'minimum v alignment': 45, 
        'minimum d alignment': 4, 
        'minimum j alignment': 30, 
        'd overlap penalty factor': 0.5
    },
    'vdj regions': {
        'VH': {
            'name': 'VH',
            'minimum identity': 0.7,
            'minimum alignment': 45
        },
        'DH': {
            'name': 'DH',
            'minimum identity': 0.7,
            'minimum alignment': 4,
            'overlap penalty factor': 0.5,
        },
        'JH': {
            'name': 'JH',
            'minimum identity': 0.8,
            'minimum alignment': 30
        }
    },
    'table': [
        {
            'collection': 'reference',
            'index': [
                { 'key': [( 'allele name', ASCENDING )], 'unique': True, 'name': 'reference allele name' },
                { 'key': [( 'region', ASCENDING )], 'unique': False, 'name': 'reference region' },
                { 'key': [( 'subgroup', ASCENDING )], 'unique': False, 'name': 'reference subgroup' },
                { 'key': [( 'gene name', ASCENDING )], 'unique': False, 'name': 'reference gene name' },
            ]
        },
        {
            'collection': 'sample',
            'index': [
                { 'key': [( 'head.id', ASCENDING ), ( 'head.library', ASCENDING )], 'unique': True, 'name': 'sample in library' },
                { 'key': [( 'head.id', ASCENDING )], 'unique': False, 'name': 'sample id' },
                { 'key': [( 'head.framed', ASCENDING )], 'unique': False, 'name': 'sample framed' },
                { 'key': [( 'head.in frame', ASCENDING )], 'unique': False, 'name': 'sample in frame' },
                { 'key': [( 'head.premature', ASCENDING )], 'unique': False, 'name': 'sample premature' },
                { 'key': [( 'head.library', ASCENDING )], 'unique': False, 'name': 'sample library' },
                { 'key': [( 'head.gapped', ASCENDING )], 'unique': False, 'name': 'sample gapped' },
                { 'key': [( 'head.valid', ASCENDING )], 'unique': False, 'name': 'sample valid' },
                { 'key': [( 'head.palindromic', ASCENDING )], 'unique': False, 'name': 'sample palindromic' },
                { 'key': [( 'head.productive', ASCENDING )], 'unique': False, 'name': 'sample productive' }
            ]
        }
    ]

}

class InvalidSampleError(Exception):
    def __init__(self, message):
        super(Exception, self).__init__(message)


class CommandLineParser(object):
    def __init__(self, node):
        self.node = node
        self.parser = ArgumentParser(**self.node['instruction'])
        self._instruction = None
        def add_argument(parser, name):
            node = self.node['prototype'][name]
            parser.add_argument(*node['flag'], **node['parameter'])

        # evaluate the type for each prototype
        for argument in self.node['prototype'].values():
            if 'type' in argument['parameter']:
                argument['parameter']['type'] = eval(argument['parameter']['type'])
                
        # add global arguments
        for argument in self.node['global']['argument']:
            add_argument(self.parser, argument)
            
        if self.sectioned:
            # Add individual command sections
            sub = self.parser.add_subparsers(**self.node['section']['instruction'])
            for action in self.node['section']['action']:
                action_parser = sub.add_parser(**action['instruction'])
                if 'argument' in action:
                    for argument in action['argument']:
                        add_argument(action_parser, argument)
                        
                # Add groups of arguments, if any.
                if 'group' in action:
                    for group in action['group']:
                        group_parser = action_parser.add_argument_group(**group['instruction'])
                        if 'argument' in group:
                            for argument in group['argument']:
                                add_argument(group_parser, argument)

    @property
    def sectioned(self):
        return 'section' in self.node and 'action' in self.node['section'] and self.node['section']['action']

    @property
    def instruction(self):
        if self._instruction == None:
            self._instruction = vars(self.parser.parse_args())
        return self._instruction

    @property
    def action(self):
        return None if 'action' not in self.instruction else self.instruction['action']

    def help(self):
        self.parser.print_help()

    @property
    def query(self):
        query = {}
        for k in [
            'library',
            'id',
            'region'
        ]:
            if k in self.instruction:
                v = self.instruction[k]
                if v is not None: query[k] = v
        for k in [
            'gapped',
            'valid',
            'framed',
            'in frame',
            'premature',
            'productive',
        ]:
            if k in self.instruction:
                if self.instruction[k] is not None:
                    if v == 'Y': query[k] = True
                    elif v == 'N': query[k] = False
        return query


class Sequence(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Sequence')
        self.pipeline = pipeline 
        self.node = node
        self._reversed = None
        self.reset()

    def __str__(self):
        buffer = StringIO()
        if self.read_frame > 0:
            buffer.write(self.nucleotide[0:self.read_frame])
            buffer.write(' ')
        for index in range(self.read_frame, self.length, 3):
            buffer.write(self.nucleotide[index:index + 3])
            buffer.write(' ')
        buffer.write('\n')
        buffer.seek(0)
        return buffer.read()

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def document(self):
        document = {}
        for k,v in self.node.items():
            if k not in [
                'codon',
                'phred',
            ]:
                document[k] = v
        return document

    def reset(self):
        if self.node is None:
            self.node = {
                'read frame': 0,
                'strand': True
            }
        for k in ['codon', 'phred']:
            if k in self.node: del self.node[k]
        self._reversed = None

    def clone(self):
        clone = {
            'nucleotide': self.nucleotide, 
            'read frame': self.read_frame,
            'strand': self.strand
        }
        if 'quality' in self.node:
            clone['quality'] = self.quality
        return Sequence(self.pipeline, clone)

    def crop(self, start, end):
        cropped = {
            'nucleotide': self.nucleotide[start:end], 
            'read frame': (3 - (start - self.read_frame)%3)%3,
            'strand': self.strand
        }
        if 'quality' in self.node:
            cropped['quality'] = self.quality[start:end]
        return Sequence(self.pipeline, cropped)

    @property
    def read_frame(self):
        return self.node['read frame']

    @read_frame.setter
    def read_frame(self, value):
        self.node['read frame'] = value
        self.reset()

    @property
    def strand(self):
        return self.node['strand']

    @strand.setter
    def strand(self, value):
        self.node['strand'] = value
        self.reset()

    @property
    def nucleotide(self):
        if 'nucleotide' in self.node:
            return self.node['nucleotide']
        else:
            return None

    @nucleotide.setter
    def nucleotide(self, value):
        self.node['nucleotide'] = value
        if not self.node['nucleotide']:
            del self.node['nucleotide']
        self.reset()

    @property
    def quality(self):
        if 'quality' in self.node:
            return self.node['quality']
        else:
            return None

    @quality.setter
    def quality(self, value):
        self.node['quality'] = value
        if not self.node['quality']:
            del self.node['quality']
        self.reset()

    @property
    def codon(self):
        if 'codon' not in self.node and self.nucleotide:
            start = self.read_frame
            end = start + 3
            codon = []
            while(not end > self.length):
                acid = self.nucleotide[start:end]
                try:
                    codon.append(self.configuration['nucleic to amino'][acid])
                except KeyError as e:
                    self.log.debug('could not resolve %s triplet to an amino acid', e)
                    if acid in self.configuration['stop codon repertoire']:
                        self.log.debug('%s is potentially a stop codon', acid)
                        codon.append('*')
                    else:
                        codon.append('X')
                start = end
                end = start + 3
            if codon:
                self.node['codon'] = ''.join(codon)
        if 'codon' in self.node:
            return self.node['codon']
        else:
            return None

    @property
    def phred(self):
        if 'phred' not in self.node and self.quality:
            self.node['phred'] = [ (ord(c) - 33) for c in self.quality ]
        return self.node['phred']

    def codon_at_frame(self, read_frame):
        codon = []
        start = read_frame
        end = start + 3
        while(not end > self.length):
            acid = self.nucleotide[start:end]
            try:
                codon.append(self.configuration['nucleic to amino'][acid])
            except KeyError as e:
                self.log.debug('could not resolve %s triplet to an amino acid', e)
                if acid in self.configuration['stop codon repertoire']:
                    self.log.debug('%s is potentially a stop codon', acid)
                    codon.append('*')
                else:
                    codon.append('X')
            start = end
            end = start + 3
        return ''.join(codon)

    @property
    def length(self):
        if self.nucleotide is not None:
            return len(self.nucleotide)
        else:
            return None

    @property
    def reversed(self):
        if self._reversed is None and self.nucleotide:
            reversed = {
                'nucleotide': ''.join([ self.configuration['complement'][b] for b in list(self.nucleotide) ][::-1]), 
                'read frame': (self.length - self.read_frame) % 3,
                'strand': not self.strand
            }
            if 'quality' in self.node:
                reversed['quality'] = self.quality[::-1]
                
            self._reversed = Sequence(self.pipeline, reversed)
        return self._reversed


class Reference(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Reference')
        self.pipeline = pipeline
        self.node = node
        
        if self.node is None:
            self.node = { 'head': {}, 'body': {} }
        
        if 'sequence' in self.body and isinstance(self.body['sequence'], dict):
            self.body['sequence'] = Sequence(self.pipeline, self.body['sequence'])
        else:
            self.body['sequence'] = Sequence(self.pipeline)

    @property
    def head(self):
        return self.node['head']

    @property
    def body(self):
        return self.node['body']

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def id(self):
        return self.head['id']

    @property
    def sequence(self):
        return self.body['sequence']

    @property
    def region(self):
        return self.head['region']

    @property
    def strain(self):
        if 'strain' in self.head:
            return self.head['strain']
        else:
            return None

    @property
    def strand(self):
        return self.body['strand']

    @property
    def framed(self):
        return self.head['framed']

    @property
    def functionality(self):
        return self.head['functionality']


class Sample(object):
    def __init__(self, pipeline, node=None, id=None, library=None):
        self.log = logging.getLogger('Sample')
        self.pipeline = pipeline
        self.node = node
        
        def load_region(node):
            if 'query' in node and isinstance(node['query'], dict):
                node['query'] = Sequence(self.pipeline, node['query'])
                
            if 'subject' in node and isinstance(node['subject'], dict):
                node['subject'] = Sequence(self.pipeline, node['subject'])

        if self.node is None:
            self.node = {
                'head': {
                    'id': id,
                    'library': library,
                    'valid': True,
                    'framed': False,
                    'gapped': False,
                    'premature': False,
                    'in frame': False,
                    'productive': False,
                    'region' : []
                }
            }
        if not self.id:
            raise InvalidSampleError('sample must have an id')
        
        if not self.key:
            raise InvalidSampleError('sample must have a valid key')
        
        if not self.library:
            raise InvalidSampleError('sample must have a library')
        
        if 'sequence' in self.node and isinstance(self.node['sequence'], dict):
            self.node['sequence'] = Sequence(self.pipeline, self.node['sequence'])
        else:
            self.node['sequence'] = Sequence(self.pipeline)
            
        if 'hit' in self.node:
            for hit in self.node['hit']:
                load_region(hit)
        else:
            self.node['hit'] = []
            
        if 'region' in self.node:
            for value in self.node['region'].values():
                load_region(value)
        else:
            self.node['region'] = {}

    def __str__(self):
        return self.diagram('default')

    @property
    def key(self):
        if 'id sha1' not in self.head:
            self.head['id sha1'] = hashlib.sha1(self.id.encode('utf8')).hexdigest()
        return self.head['id sha1']

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def document(self):
        def transform_to_document(node):
            if isinstance(node, list):
                return [ transform_to_document(v) for v in node ]
                
            elif isinstance(node, dict):
                for key, value in node.items():
                    node[key] = transform_to_document(value)
                return node
                
            elif isinstance(node, Sequence):
                return node.document
                
            else: return node
            
        document = transform_to_document(self.node)
        if document['hit']:
            document['matched'] = True
        else:
            document['matched'] = False
            del document['hit']
            del document['region']
        return document

    @property
    def head(self):
        return self.node['head']

    @property
    def id(self):
        return self.head['id']

    @id.setter
    def id(self, value):
        self.head['id'] = value

    @property
    def length(self):
        return self.sequence.length

    @property
    def strand(self):
        return self.node['strand']

    @strand.setter
    def strand(self, value):
        self.node['strand'] = value

    @property
    def library(self):
        return self.head['library']

    @library.setter
    def library(self, value):
        self.head['library'] = value

    @property
    def sequence(self):
        return self.node['sequence']

    @property
    def valid(self):
        return self.head['valid']

    @valid.setter
    def valid(self, value):
        self.head['valid'] = value

    @property
    def hit(self):
        return self.node['hit']

    @property
    def region(self):
        if 'region' not in self.node:
            self.node['region'] = {}
        return self.node['region']

    @property
    def matched(self):
        return len(self.hit) > 0

    @property
    def framed(self):
        return self.head['framed']

    @property
    def gapped(self):
        return self.head['gapped']

    @property
    def in_frame(self):
        return self.head['in frame']

    @property
    def premature(self):
        return self.head['premature']

    @property
    def productive(self):
        return self.head['productive']

    @property
    def fasta(self):
        buffer = []
        buffer.append('>' + self.id)
        begin = 0
        end = 0
        while end < self.sequence.length:
            end = min(begin + self.configuration['constant']['fasta line length'], self.sequence.length)
            buffer.append(self.sequence.nucleotide[begin:end])
            begin = end
        return '\n'.join(buffer)

    @property
    def fastq(self):
        return '{}\n{}\n+\n{}'.format(self.id, self.sequence.nucleotide, self.sequence.quality)

    def add_hit(self, hit):
        if hit is not None:
            hit['query start'] -= 1
            if 'uuid' not in hit: hit['uuid'] = str(uuid.uuid4())
            self.hit.append(hit)

    def reverse(self):
        self.node['sequence'] = self.node['sequence'].reversed

    def expand(self):
        if 'hit' in self.node:
            for hit in self.node['hit']:
                if 'compressed hit' in hit:
                    match = self.configuration['expression']['expand hit'].search(hit['compressed hit'])
                    if match:
                        parsed = dict(((k.replace('_', ' '),v) for k,v in match.groupdict().items()))
                        for k,v in parsed.items():
                            if k in [
                                'query start',
                                'query end',
                                'subject start',
                                'subject end',
                                'gap openings',
                                'gaps',
                                'mismatch',
                                'alignment length'
                            ]:
                                try:
                                    hit[k] = int(parsed[k])
                                except ValueError as e:
                                    self.log.warning('could not parse value %s for %s as int', parsed[k], k)
                            elif k in [
                                'identical',
                                'bit score',
                                'evalue',
                            ]:
                                try:
                                    hit[k] = float(parsed[k])
                                except ValueError:
                                    self.log.warning('could not parse value %s for %s as float', parsed[k], k)
                            else:
                                hit[k] = v

    def compress(self):
        if 'hit' in self.node:
            for hit in self.node['hit']:
                try:
                    hit['compressed hit'] = self.configuration['expression']['igblast compressed hit'].format(**hit)
                except KeyError as e:
                    self.log.error('could not compress hit because %s was missing', e)

    def analyze(self, strain):
        self.reset()
        if self.valid: self._load_reference()
        if self.valid: self._pick_jh_region(strain)
        if self.valid: self._pick_vh_region(strain)
        if self.valid: self._check_for_aligned_frames()
        if self.valid: self._pick_dh_region()
        if self.valid: self._identify_v_j_junction()
        if self.valid: self._identify_v_d_junction()
        if self.valid: self._identify_d_j_junction()
        if self.valid: self._check_for_stop_codon()
        if self.valid: self._identify_cdr3()
        if self.valid: self._check_productive()
        if self.valid: self._identify_chewback()
        if self.valid: self._index_regions()
        if self.valid: self._analyze_quality()

    def to_json(self):
        def handler(o):
            result = None
            if isinstance(o, datetime):
                result = o.isoformat()
            if isinstance(o, ObjectId):
                result = str(o)
            if isinstance(o, set):
                result = list(o)
            if isinstance(o, Sequence):
                result = o.node
            return result
            
        return json.dumps(self.document, sort_keys=True, ensure_ascii=False, indent=4, default=handler)

    def to_summary(self):
        if self.hit:
            buffer = []
            buffer.append(self.id)
            for r in ('VH', 'DH', 'JH'):
                if r in self.region:
                    buffer.append('{} : {}'.format(r, self.region[r]['subject id']))
            
            if self.productive:
                buffer.append('productive')
            else:
                if self.framed:
                    buffer.append('in frame' if self.in_frame else 'out of frame')
                    
                if self.premature:
                    buffer.append('premature')
                
            buffer.append('Q : {:.4}'.format(self.head['average phred']))
            return ' | '.join(buffer)
        else:
            return ''

    def view(self, profile):
        diagram = Diagram(self.pipeline, self, profile)
        print(diagram.draw())

    def info(self, profile):
        print(self.to_json())

    def invalidate(self, message):
        self.valid = False
        self.node['error'] = message
        self.log.warning('sample %s invalidated because %s', self.id, message)

    def invalidate_hit(self, hit, message):
        hit['valid'] = False
        hit['error'] = message
        self.log.warning('%dbp %s hit to %s on %s invalidated because %s', 
        hit['alignment length'], hit['region'], hit['subject id'], self.id, message)

    def _set_frame(self, hit):
        if hit['framed']:
            self.head['framed'] = True
            self.node['framed by'] = hit['subject id']
            self.sequence.read_frame = (hit['query start'] + hit['subject'].read_frame) % 3
            self.strand = hit['subject strand']
            self._assign_frame()
            self.log.debug('frame set for %s by %s', self.id, hit['subject id'])

    def reset(self):
        self.head['valid'] = True
        self.head['gapped'] = False
        self.head['framed'] = False
        self.head['premature'] = False
        self.head['in frame'] = False
        self.head['productive'] = False
        self.head['palindromic'] = False
        self.head['region'] = []
        for hit in self.hit:
            hit['valid'] = True
            hit['in frame'] = False
            hit['picked'] = False
            hit['framed'] = False
            hit['gapped'] = False
            hit['score'] = hit['bit score']
            if 'gap openings' in hit and hit['gap openings'] > 0:
                hit['gapped'] = True
                self.head['gapped'] = True
                self.invalidate_hit(hit, 'hit contains a gapped alignment')
                
            if hit['valid']:
                hit['query'] = self.sequence.crop(hit['query start'], hit['query end'])

    def _load_reference(self):
        for hit in self.hit:
            if hit['valid']:
                if 'subject id' in hit:
                    reference = self.pipeline.reference_for(hit['subject id'])
                    if reference:
                        if hit['subject strand'] == reference.sequence.strand:
                            ref = reference.sequence
                            hit['subject start'] -= 1
                        else:
                            ref = reference.sequence.reversed
                            hit['subject start'] = ref.length - hit['subject start']
                            hit['subject end'] = ref.length - hit['subject end'] + 1
            
                        hit['subject'] = ref.crop(hit['subject start'], hit['subject end'])
                        hit['framed'] = reference.framed
                        hit['functionality'] = reference.functionality
                        if reference.strain: hit['strain'] = reference.strain
                            
                        # only DH regions are allowed to align to the oposite strand
                        if hit['region'] != 'DH' and hit['subject strand'] != reference.sequence.strand:
                            self.invalidate_hit(hit, 'hit aligns to the wrong strand')

    def _assign_frame(self):
        for hit in self.hit:
            if hit['valid']:
                hit['query'].read_frame = 2 - (hit['query start'] - self.sequence.read_frame - 1) % 3
                if hit['framed']:
                    if hit['query'].read_frame == hit['subject'].read_frame:
                        hit['in frame'] = True

    def _pick_region(self, region, candidate, strain=None):
        picked = False
        if region in self.configuration['vdj regions']:
            region_node = self.configuration['vdj regions'][region]
            top = None
            search = { 'valid': [] }
            for hit in candidate:
                if (hit['valid'] and 
                    hit['region'] == region_node['name'] and
                    hit['identical'] >= region_node['minimum identity'] and 
                    hit['alignment length'] >= region_node['minimum alignment']):
                    search['valid'].append(hit)
                    
            if search['valid']:
                top = search['valid']
                top.sort(reverse=True, key=lambda x: x['score'])
                top = [ hit for hit in top if hit['score'] == top[0]['score'] ]
                
                if len(top) > 1:
                    search['framed'] = []
                    sample_frame = None if not self.framed else self.sequence.read_frame
                    for hit in top:
                        frame = hit['subject'].read_frame if hit['framed'] else sample_frame
                        if frame is not None:
                            subject = hit['subject'].codon_at_frame(frame)
                            query = hit['query'].codon_at_frame(frame)
                            hit['codon mismatch'] = 0
                            for index,c in enumerate(subject):
                                if c != query[index]:
                                    hit['codon mismatch'] += 1
                            search['framed'].append(hit)
                            
                    if search['framed']:
                        top = search['framed']
                        top.sort(key=lambda x: x['codon mismatch'])
                        top = [ hit for hit in top if hit['codon mismatch'] == top[0]['codon mismatch'] ]
    
                if len(top) > 1 and strain is not None:
                    search['strained'] = [ hit for hit in top if 'strain' in hit and hit['strain'] == strain ]
                    if search['strained']:
                        top = search['strained']
                        
                if len(top) > 1:
                    search['functional'] = [ hit for hit in top if 'functionality' in hit and hit['functionality'] == 'F' ]
                    if search['functional']:
                        top = search['functional']
                        
            if top:
                self.region[region] = top[0]
                for hit in top: hit['picked'] = True
                if not self.framed: self._set_frame(hit)
                picked = True
                # if not picked: self.invalidate('no satisfactory %s region match found', region_node['name'])
        return picked

    def _pick_jh_region(self, strain):
        return self._pick_region('JH', self.hit, strain)

    def _pick_vh_region(self, strain):
        return self._pick_region('VH', self.hit, strain)

    def _pick_dh_region(self):
        region_node = self.configuration['vdj regions']['DH']
        top = None
        search = { 'valid': [], 'map': {} }
        for hit in self.hit:
            if hit['region'] == 'DH':
                search['map'][hit['uuid']] = hit
                if hit['valid']:
                    search['valid'].append(hit)
                    
        if search['valid']:
            if len(search['valid']) > 0:
                search['trimmed'] = []
                for hit in search['valid']:
                    trimmed = { 'overlap': 0 }
                    for k,v in hit.items():
                        if k not in ('query', 'subject', 'overlap'):
                            trimmed[k] = v
                            
                    reference = self.pipeline.reference_for(hit['subject id'])
                    
                    # check for overlap with VH
                    if 'VH' in self.region and self.region['VH']['query end'] > hit['query start']:
                        trimmed['query start'] = self.region['VH']['query end']
                        overlap = trimmed['query start'] - hit['query start']
                        trimmed['subject start'] += overlap
                        trimmed['overlap'] += overlap
                        
                    # check for overlap with JH
                    if 'JH' in self.region and hit['query end'] > self.region['JH']['query start']:
                        trimmed['query end'] = self.region['JH']['query start']
                        overlap = hit['query end'] - trimmed['query end']
                        trimmed['subject end'] -= overlap
                        trimmed['overlap'] += overlap
                        
                    # correct the score for the overlap and trim the sequences
                    trimmed['score'] = hit['bit score'] - trimmed['overlap'] * region_node['overlap penalty factor']
                    trimmed['alignment length'] = trimmed['query end'] - trimmed['query start']
                    trimmed['query'] = self.sequence.crop(trimmed['query start'], trimmed['query end'])
    
                    if trimmed['subject strand'] == reference.sequence.strand:
                        trimmed['subject'] = reference.sequence.crop(trimmed['subject start'], trimmed['subject end'])
                    else:
                        trimmed['subject'] = reference.sequence.reversed.crop(trimmed['subject start'], trimmed['subject end'])
                                        
                    # filter DH hits that match the criteria
                    if trimmed['query end'] > trimmed['query start']:
                        similar = 0
                        different = 0
                        longest = 0
                        stretch = 0
                        for index in range(trimmed['query'].length):
                            if trimmed['subject'].nucleotide[index] == trimmed['query'].nucleotide[index]:
                                similar += 1
                                stretch += 1
                                longest = max(longest, stretch)
                            else:
                                different += 1
                                stretch = 0
                        trimmed['identical'] = float(similar) / float(trimmed['alignment length'])
                        original = search['map'][trimmed['uuid']]
                        for k,v in trimmed.items(): original[k] = v
                        if longest >= region_node['minimum alignment']:
                            search['trimmed'].append(original)
                            
                if search['trimmed']:
                    top = search['trimmed']
        if top:
            return self._pick_region('DH', top)
        else:
            return False

    def _identify_cdr3(self):
        if 'VH' in self.region and 'JH' in self.region:
            start = None
            end = None
            
            offset = None
            # look for the most upstream cycteine on the VH region
            for index,codon in enumerate(reversed(self.region['VH']['query'].codon)):
                if codon == 'C':
                    offset = self.region['VH']['query'].read_frame + (len(self.region['VH']['query'].codon) - index - 1) * 3
                    break
                    
            if offset is not None:
                start = self.region['VH']['query start'] + offset
                
            offset = None
            # look for the most downstream tryptophan JH region
            # followed by GG, followed by either C or T
            for index,codon in enumerate(self.region['JH']['query'].codon):
                if codon == 'W':
                    o = self.region['JH']['query'].read_frame + (index + 1) * 3
                    suffix = self.region['JH']['query'].nucleotide[o:o + 3]
                    if suffix[0:2] == 'GG':
                        if suffix[2] == 'T' or suffix[2] == 'C':
                            offset = o
                            break

            if offset is not None:
                end = self.region['JH']['query start'] + offset
                
            if start and end:
                self.region['CDR3'] = {
                    'subject strand': self.strand,
                    'picked': True,
                    'subject id': 'CDR3',
                    'query start': start,
                    'query end': end,
                    'region': 'CDR3',
                    'query': self.sequence.crop(start, end)
                }

    def _check_for_stop_codon(self):
        if '*' in self.sequence.codon:
            self.head['premature'] = True
        else:
            self.head['premature'] = False
        return not self.head['premature']

    def _check_productive(self):
        if (self.valid and
            not self.premature and
            self.in_frame and 
            'VH' in self.region and
            'JH' in self.region):
            self.head['productive'] = True

    def _identify_chewback(self):
        for hit in self.hit:
            if hit['valid']:
                if 'subject id' in hit:
                    reference = self.pipeline.reference_for(hit['subject id'])
                    if reference:
                        if hit['subject strand'] == reference.sequence.strand:
                            ref = reference.sequence
                        else:
                            ref = reference.sequence.reversed
            
                        if hit['region'] == 'VH' or hit['region'] == 'DH':
                            if hit['subject end'] < ref.length:
                                hit['3 chew'] = ref.crop(hit['subject end'], ref.length)
                                
                        if hit['region'] == 'JH' or hit['region'] == 'DH':
                            if hit['subject start'] > 0:
                                hit['5 chew'] = ref.crop(0, hit['subject start'])

    def _analyze_quality(self):
        for name, region in self.region.items():
            self.head['{} length'.format(name)] = region['query'].length
            region['average phread'] = float(sum(region['query'].phred)) / float((region['query'].length))
        self.head['average phred'] = float(sum(self.sequence.phred)) / float((self.sequence.length))

    def _check_for_aligned_frames(self):
        if 'JH' in self.region and 'VH' in self.region:
            if self.region['VH']['in frame'] and self.region['JH']['in frame']:
                self.head['in frame'] = True

    def _index_regions(self):
        self.head['region'] = []
        for hit in self.hit:
            if hit['valid'] and hit['picked']:
                self.head['region'].append(hit['subject id'])
        for region in [
            'CDR3',
            'V-D',
            'D-J',
            'V-J'
        ]: 
            if region in self.region:
                self.head['region'].append(region)
        self.head['region'] = list(set(self.head['region']))

    def _check_palindrome(self, junction, left, right):
        if junction['query'].length > 0:
            junction['palindrome'] = [ 'N' ] * junction['query'].length
            i = 0
            while i < left.length and i < junction['query'].length:
                mirror = -(i + 1)
                nucleotide = junction['query'].nucleotide[i]
                palindrome = left.nucleotide[mirror]
                if self.configuration['complement'][nucleotide] == palindrome:
                    junction['palindrome'][i] = 'P'
                    i += 1
                else:
                    break
            i = 0
            while i < right.length and i < junction['query'].length:
                mirror = -(i + 1)
                nucleotide = right.nucleotide[i]
                palindrome = junction['query'].nucleotide[mirror]
                if self.configuration['complement'][nucleotide] == palindrome:
                    junction['palindrome'][mirror] = 'P'
                    i += 1
                else:
                    break
            junction['palindrome'] = ''.join(junction['palindrome'])
            p = junction['palindrome'].count('P')
            n = junction['palindrome'].count('N')
            if p > 0: self.head['palindromic'] = True
            if n > 0:
                if p > 0:
                    junction['palindrome ratio'] = float(n) / float(p)
                else:
                    junction['palindrome ratio'] = 1.0
            else:
                junction['palindrome ratio'] = 0.0

    def _identify_v_d_junction(self):
        if 'DH' in self.region and 'VH' in self.region:
            if self.region['VH']['query end'] < self.region['DH']['query start']:
                start = self.region['VH']['query end']
                end = self.region['DH']['query start']
                junction = {
                    'subject strand': self.strand,
                    'picked': True,
                    'subject id': 'V-D',
                    'query start': start,
                    'query end': end,
                    'region': 'V-D',
                    'query': self.sequence.crop(start, end)
                }
                self._check_palindrome(junction, self.region['VH']['query'], self.region['DH']['query'])
                self.region['V-D'] = junction    
        return True

    def _identify_d_j_junction(self):
        if 'DH' in self.region and 'JH' in self.region:
            if self.region['JH']['query start'] > self.region['DH']['query end']:
                start = self.region['DH']['query end']
                end = self.region['JH']['query start']
                junction = {
                    'subject strand': self.strand,
                    'picked': True,
                    'subject id': 'D-J',
                    'query start': start,
                    'query end': end,
                    'region': 'D-J',
                    'query': self.sequence.crop(start, end)
                }
                self._check_palindrome(junction, self.region['DH']['query'], self.region['JH']['query'])
                self.region['D-J'] = junction
        return True

    def _identify_v_j_junction(self):
        if 'DH' not in self.region and 'JH' in self.region and 'VH' in self.region:
            if self.region['JH']['query start'] > self.region['VH']['query end']:
                start = self.region['VH']['query end']
                end = self.region['JH']['query start']
                junction = {
                    'subject strand': self.strand,
                    'picked': True,
                    'subject id': 'V-J',
                    'query start': start,
                    'query end': end,
                    'region': 'V-J',
                    'query': self.sequence.crop(start, end)
                }
                self._check_palindrome(junction, self.region['VH']['query'], self.region['JH']['query'])
                self.region['V-J'] = junction
        return True


class Diagram(object):
    def __init__(self, pipeline, sample, profile):
        self.log = logging.getLogger('Diagram')
        self.pipeline = pipeline
        self.sample = sample
        self.node = {
            'track offset': {},
            'width': {},
            'query': {},
            'track': [],
            'pattern': {
                'width': 0,
                'title': [],
                'phrase': [],
                'feature': [],
            },
        }

        if profile in self.configuration['profile'] and 'diagram' in self.configuration['profile'][profile]:
            p = self.configuration['profile'][profile]['diagram']
        else:
            p = self.configuration['profile']['default']['diagram']
            
        for k,v in p['track'].items():
            self.query[k] = v
        
        for track in sample.hit:
            if all([k in track and track[k] == v for k,v in self.query.items()]):
                self.track.append(track)

        for track in sample.region.values():
            if track['subject id'] in set(['CDR3', 'V-J', 'V-D', 'D-J']):
                if 'uuid' not in track:
                    track['uuid'] = str(uuid.uuid4())
                self.track.append(track)
            
        for k in p['feature']:
            if k in self.configuration['diagram']['prototype']:
                feature = dict(self.configuration['diagram']['prototype'][k])
                if feature['width'] == 'auto':
                    feature['width'] = 0
                    for track in self.track:
                        if feature['value'] in track:
                            feature['width'] = max(feature['width'], len(track[feature['value']]))
                self.pattern['feature'].append(feature)
                self.pattern['phrase'].append('{{: <{:}}}'.format(feature['width']))
                self.pattern['title'].append(feature['title'])
                self.pattern['width'] += max(feature['width'], len(feature['title']))
        self.pattern['width'] += (len(self.pattern['phrase']) - 1)
        self.pattern['phrase'] = ' '.join(self.pattern['phrase'])
        self.width['sample read frame'] = self.sample.sequence.read_frame
        self.width['table'] = self.pattern['width']
        self.width['table padding'] = 2
        self.width['diagram start'] = self.width['table'] + self.width['table padding']
        
        
        for track in self.track:
            self.node['track offset'][track['uuid']] = self._find_track_offset(track)
                
        

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def track(self):
        return self.node['track']

    @property
    def pattern(self):
        return self.node['pattern']

    @property
    def query(self):
        return self.node['query']

    @property
    def framed(self):
        return self.sample.framed

    @property
    def gap(self):
        return ' ' if self.sample.framed else ''

    @property
    def width(self):
        return self.node['width']

    def format_track_features(self, track):
        values = []
        for feature in self.pattern['feature']:
            if feature['value'] in track:
                values.append(feature['format'](track[feature['value']]))
            else:
                values.append('')
        return self.pattern['phrase'].format(*values)

    def _draw_summary(self, buffer):
        b = []
        b.append(self.sample.library)
        b.append(self.sample.id)
        for r in ('VH', 'DH', 'JH'):
            if r in self.sample.region:
                b.append('{} : {}'.format(r, self.sample.region[r]['subject id']))
        
        if self.sample.productive:
            b.append('productive')
        else:
            if self.sample.framed:
                b.append('in frame' if self.sample.in_frame else 'out of frame')
                
            if self.sample.premature:
                b.append('premature')
        b.append('Q {:.4}'.format(self.sample.head['average phred']))
        
        buffer.write(' | '.join(b))
        buffer.write('\n')

    def _draw_track_title(self, buffer):
        buffer.write(self.pattern['phrase'].format(*(self.pattern['title'])))

    def _draw_coordinates(self, buffer):
        buffer.write(' ' * self.width['diagram start'])
        if self.width['sample read frame'] > 0:
            buffer.write('{: <{}}'.format(0, self.width['sample read frame']))
            buffer.write(self.gap)
            
        for index in range(self.width['sample read frame'], self.sample.sequence.length, 3):
            buffer.write('{: <3}'.format(index))
            buffer.write(self.gap)
        buffer.write('\n')

    def _draw_nucleotide_sequence(self, buffer):
        buffer.write(' ' * self.width['diagram start'])
        if self.width['sample read frame'] > 0:
            buffer.write(self.sample.sequence.nucleotide[0:self.width['sample read frame']])
            buffer.write(self.gap)
        for index in range(self.width['sample read frame'], self.sample.sequence.length, 3):
            buffer.write(self.sample.sequence.nucleotide[index:index + 3])
            buffer.write(self.gap)
        buffer.write('\n')

    def _draw_quality_sequence(self, buffer):
        buffer.write(' ' * self.width['diagram start'])
        if self.width['sample read frame'] > 0:
            buffer.write(self.sample.sequence.quality[0:self.width['sample read frame']])
            buffer.write(self.gap)
        for index in range(self.width['sample read frame'], self.sample.sequence.length, 3):
            buffer.write(self.sample.sequence.quality[index:index + 3])
            buffer.write(self.gap)
        buffer.write('\n')

    def _draw_codon_sequence(self, buffer):
        if self.sample.framed:
            buffer.write(' ' * self.width['diagram start'])
            if self.width['sample read frame'] > 0:
                buffer.write(' ' * self.width['sample read frame'])
                buffer.write(self.gap)
            for codon in self.sample.sequence.codon:
                buffer.write('{: <3}'.format(codon))
                buffer.write(self.gap)
            buffer.write('\n')

    def _draw_track_without_subject(self, buffer, track):
        offset = self.node['track offset'][track['uuid']]
        if 'subject' not in track and 'query' in track:
            if offset > 0:
                buffer.write(' ' * offset)
            if track['query'].read_frame > 0:
                buffer.write('-' * track['query'].read_frame)
                # buffer.write(track['query'].nucleotide[0:track['query'].read_frame])
                buffer.write(self.gap)
                
            for index in range(track['query'].read_frame, track['query'].length, 3):
                buffer.write('-' * len(track['query'].nucleotide[index:index + 3]))
                # buffer.write(track['query'].nucleotide[index:index + 3])
                buffer.write(self.gap)
            buffer.write('\n')
            
            if 'palindrome' in track and 'P' in track['palindrome']:
                buffer.write(' ' * self.width['diagram start'])
                if offset > 0: buffer.write(' ' * offset)
                if track['query'].read_frame > 0:
                    buffer.write(track['palindrome'][0:track['query'].read_frame])
                    buffer.write(self.gap)
                    
                for index in range(track['query'].read_frame, track['query'].length, 3):
                    buffer.write(track['palindrome'][index:index + 3])
                    buffer.write(self.gap)
                buffer.write('\n')
                
            if False and self.sample.framed and track['query'].codon:
                buffer.write(' ' * self.width['diagram start'])
                if offset > 0: buffer.write(' ' * offset)
                    
                if track['query'].read_frame > 0:
                    buffer.write(' ' * track['query'].read_frame)
                    buffer.write(self.gap)
                    
                for codon in track['query'].codon:
                    buffer.write('{: <3}'.format(codon))
                    buffer.write(self.gap)
                buffer.write('\n')

    def _draw_track_with_subject(self, buffer, track):
        offset = self.node['track offset'][track['uuid']]
        if 'subject' in track and 'query' in track:
            subject = track['subject'].clone()
            subject.read_frame = track['query'].read_frame
            
            mask = []
            for index,n in enumerate(subject.nucleotide):
                mask.append('-' if n == track['query'].nucleotide[index] else n)
            mask = ''.join(mask)
                
            if offset > 0: buffer.write(' ' * offset)
            if subject.read_frame > 0:
                buffer.write(mask[0:subject.read_frame])
                buffer.write(self.gap)
                
            for index in range(subject.read_frame, subject.length, 3):
                buffer.write(mask[index:index + 3])
                buffer.write(self.gap)
            buffer.write('\n')
            
            if self.sample.framed and subject.codon:
                display = False
                mask = []
                for index,c in enumerate(subject.codon):
                    if c == track['query'].codon[index]: mask.append('-')
                    else:
                        mask.append(c)
                        display = True
                if display:
                    buffer.write(' ' * self.width['diagram start'])
                    if offset > 0: buffer.write(' ' * offset)
                        
                    if subject.read_frame > 0:
                        buffer.write(' ' * subject.read_frame)
                        buffer.write(self.gap)
                        
                    for codon in mask:
                        buffer.write('{: <3}'.format(codon))
                        buffer.write(self.gap)
                    buffer.write('\n')

    def _find_track_offset(self, track):
        result = 0
        if track['query start'] > 0:
            result = track['query start']
            if self.sample.framed:
                if self.width['sample read frame'] > 0:
                    result += int((track['query start'] - self.width['sample read frame']) / 3) + 1
                else:
                    result += (int(track['query start'] / 3))
        return result

    def _draw_track(self, buffer, track):
        buffer.write(self.format_track_features(track))
        buffer.write(' ' * self.width['table padding'])
        self._draw_track_without_subject(buffer, track)
        self._draw_track_with_subject(buffer, track)

    def draw(self):
        buffer = StringIO()
        buffer.write('\n')
        self._draw_track_title(buffer)
        buffer.write(' ' * self.width['table padding'])
        self._draw_summary(buffer)
        self._draw_coordinates(buffer)
        self._draw_nucleotide_sequence(buffer)
        self._draw_quality_sequence(buffer)
        self._draw_codon_sequence(buffer)
        for track in self.track:
            self._draw_track(buffer, track)
        buffer.seek(0)
        return buffer.read()


class Block(object):
    def __init__(self, pipeline):
        self.log = logging.getLogger('Block')
        self.pipeline = pipeline
        self.reset()

    def __str__(self):
        buffer = StringIO()
        for sample in self.buffer:
            if sample.valid:
                buffer.write(str(sample))
                buffer.write('\n')
        buffer.seek(0)
        return buffer.read()

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def buffer(self):
        return self.node['buffer']

    @property
    def document(self):
        return [ sample.document for sample in self.buffer ]

    @property
    def lookup(self):
        return self.node['lookup']

    @property
    def size(self):
        return len(self.buffer)

    @property
    def empty(self):
        return self.size == 0

    @property
    def full(self):
        return not self.size < self.configuration['constant']['buffer size']

    def reset(self):
        self.node = {
            'buffer': [],
            'lookup': {}
        }

    def add(self, sample):
        #if sample is not None and sample.sequence.length > 0:
        if sample.id not in self.lookup:
            self.buffer.append(sample)
            self.lookup[sample.id] = sample
        else:
            self.log.error('%s already present', sample.id)

    def fill(self, library, strain):
        gapped = False
        if self.read(library):
            buffer = self.search()
            if buffer:
                self.parse_igblast(buffer)
                for sample in self.buffer:
                    sample.analyze(strain)
                    if sample.gapped:
                        gapped = True
            if gapped:
                buffer = self.search_gapped()
                if buffer:
                    self.parse_igblast_gapped(buffer)
        return not self.empty

    def read(self, library):
        self.reset()
        state = 0
        sample = None
        for line in sys.stdin:
            if not line:
                break
            else:
                line = line.strip()
                if line:
                    if state == 0 and line[0] == '@':
                        try:
                            sample = Sample(self.pipeline, id=line, library=library)
                        except InvalidSampleError as e:
                            self.log.warning(e)
                        state = 1
                            
                    elif state == 1:
                        if sample is not None:
                            sample.sequence.nucleotide = line
                        state = 2
                        
                    elif state == 2:
                        state = 3
                            
                    elif state == 3:
                        if sample is not None:
                            sample.sequence.quality = line
                            self.add(sample)
                            sample = None
                            state = 0
                            if self.full:
                                break
        return not self.empty

    def search(self):
        result = None
        command = self.configuration['command']['igblast']
        process = Popen(
            args=command['arguments'],
            cwd=command['cwd'],
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate(input=self.to_fasta().read().encode('utf8'))
        if output: result = StringIO(output.decode('utf8'))
        return result

    def search_gapped(self):
        result = None
        command = self.configuration['command']['igblast.gapped']
        process = Popen(
            args=command['arguments'],
            cwd=command['cwd'],
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        fasta = self.to_fasta({'gapped': True})
        fasta.seek(0)
        output, error = process.communicate(input=fasta.read().encode('utf8'))
        if output: result = StringIO(output.decode('utf8'))
        return result

    def to_fasta(self, query=None):
        buffer = StringIO()
        for sample in self.buffer:
            if query is None or all([k in sample.head and sample.head[k] == v for k,v in query.items()]):
                buffer.write('>')
                buffer.write(sample.id)
                buffer.write('\n')
                begin = 0
                end = 0
                while end < sample.sequence.length:
                    end = min(begin + self.configuration['constant']['fasta line length'], sample.sequence.length)
                    buffer.write(sample.sequence.nucleotide[begin:end])
                    buffer.write('\n')
                    begin = end
        buffer.seek(0)
        return buffer

    def parse_igblast_gapped(self, buffer):
        buffer.seek(0)
        results = []
        state = 0
        for line in buffer:
            if state == 0 and line[0:6] == 'Query=':
                record = {
                    'id': line.strip()[6:].strip() ,
                    'sequence': '',
                }
                record['sample'] = self.lookup[record['id']]
                state = 1
                
            elif state == 1 and line[0:10] == 'Alignments':
                state = 2
            elif state == 2:
                match = configuration['expression']['gapped sequence'].search(line)
                if match:
                    g = match.groupdict()
                    if 'offset' not in record:
                        record['offset'] = int(g['offset']) - 1
                    record['sequence'] += g['sequence']
                else:
                    if line[0:6] == 'Lambda':
                        state = 0
                        results.append(record)
                        
        for record in results:
            one = (' ' * record['offset']) + record['sequence']
            two = record['sample'].sequence.nucleotide
            o = 0
            t = 0
            r = []
            while o < len(one) or t < len(two):
                if not(o < len(one)):
                  r.append(two[t])
                  o += 1
                  t += 1
                elif not(o < len(two)):
                  r.append(one[o])
                  o += 1
                  t += 1
                elif one[o] == ' ':
                  r.append(two[t])
                  o += 1
                  t += 1
                elif one[o] == '-':
                  r.append('-')
                  o += 1
                elif two[t] == '-':
                  r.append('-')
                  t += 1
                else:
                    if one[o] == two[t]:
                        r.append(one[o])
                        o += 1
                        t += 1
                    else:
                        r.append('*')
                        o += 1
                        t += 1
            record['sample'].head['gapped sequence'] = ''.join(r)

    def parse_igblast(self, buffer):
        def parse_igblast_hit(line):
            hit = None
            match = self.configuration['expression']['igblast hit'].search(line)
            if match:
                hit = dict(((k.replace('_', ' '),v) for k,v in match.groupdict().items()))
                if hit['subject strand'] == 'plus':
                    hit['subject strand'] = True
                elif hit['subject strand'] == 'minus':
                    hit['subject strand'] = False
                else:
                    self.log.warning('could not parse value %s for as a strand setting to plus', hit['subject strand'])
                    hit['subject strand'] = True
                        
                for k in [
                    'identical',
                    'bit score',
                    'evalue',
                ]: 
                    try: hit[k] = float(hit[k])
                    except ValueError:
                        self.log.warning('could not parse value %s for %s as int', hit[k], k)
                        hit[k] = None
                    
                for k in [
                    'subject end',
                    'gap openings',
                    'subject start',
                    'mismatch',
                    'query start',
                    'alignment length',
                    'gaps',
                    'query end',
                ]:
                    try: hit[k] = int(hit[k])
                    except ValueError as e:
                        self.log.warning('could not parse value %s for %s as int', hit[k], k)
                        hit[k] = None

                # fix the region for heavy chain
                if 'region' in hit: hit['region'] += 'H'
                
                # switch from % to fraction
                if 'identical' in hit: hit['identical'] /= 100.0
            return hit

        state = 0
        sample = None
        buffer.seek(0)
        for line in buffer:
            line = line.strip()
            
            if line.startswith('# IGBLASTN '):
                # this is the start of a new record
                if state == 2:
                    sample.invalidate('no alignment results found')
                state = 1
                sample = None
                
            elif state == 1 and line.startswith('# Query: '):
                # this is the query id for the record
                id = line[9:]
                if id in self.lookup:
                    sample = self.lookup[id]
                    state = 2
                else:
                    self.log.error('could not locate %s in buffer', id)
                    state = 0
                    
            elif state == 2 and line.startswith(self.configuration['expression']['igblast reversed query']):
                # this means the hits will be for the reverse complement strand
                sample.reverse()
                
            elif state == 2 and line.startswith('# Hit table '):
                # this means the hit table is on the next line
                state = 3
                
            elif state == 3 and not line.startswith('#'):
                hit = parse_igblast_hit(line)
                if hit is not None:
                    sample.add_hit(hit)

    def view(self, profile):
        for sample in self.buffer:
            sample.view(profile)

    def info(self, prpfile):
        for sample in self.buffer:
            sample.info(profile)

    def simulate(self, json, alignment, profile):
        for sample in self.buffer:
            if json: sample.info(profile)
            if alignment: sample.view(profile)


class Pipeline(object):
    def __init__(self):
        self.log = logging.getLogger('Pipeline')
        self.configuration = configuration
        self._connection = None
        self._reference_sequence = {}
        self._accession = {}
        self._stop_codon_feature = None
        self._strains = None
        self._load_configuration()

    def _load_configuration(self):
        # load a reverse complement table for ambiguity code
        self.configuration['complement'] = {}
        for k,v in self.configuration['iupac nucleic acid notation'].items():
            self.configuration['complement'][k] = v['reverse']
        
        # load collection of codons possibly encoding a stop codon
        stop = [ k for k,v in self.configuration['nucleic to amino'].items() if v == '*' ]
        motif = self.motif_for(stop)
        self.configuration['stop codon repertoire'] = motif['space']

    @property
    def connection(self):
        if self._connection is None:
            try:
                self._connection = MongoClient('mongodb://localhost/somatic')
                #self._connection = MongoClient('mongodb://somatic:fPWZq8nCVizzHloVDhs=@albireo.bio.nyu.edu/somatic')
            except pymongo.errors.ConnectionFailure as e:
                self.log.error('failed to establish connection %s', e)
            else:
                self.log.debug('connection established')
        return self._connection

    @property
    def database(self):
        if self.connection is not None:
            return self.connection['somatic']
        else:
            return None

    @property
    def strains(self):
        if self._strains == None:
            self._strains = list(self.database['reference'].distinct('strain'))
        return self._strains

    def close(self):
        if self._connection is not None:
            self._connection.close()

    def motif_for(self, codons):
        def expand(motif, space=[ '' ]):
            position = len(space[0])
            if position < len(motif):
                next = []
                for option in motif[position]:
                    for element in space:
                        next.append(element + option)
                return expand(motif, next)
            else:
                return space
    
        feature = { 'codon': {}, 'space': set() }
        for c in codons:
            feature['codon'][c] = { 'motif':[] }
        for triplet,codon in feature['codon'].items():
            for nucleotide in triplet:
                motif = []
                codon['motif'].append(motif)
                for k,n in configuration['iupac nucleic acid notation'].items():
                    if nucleotide in n['option']:
                        motif.append(k)
            codon['possible'] = expand(codon['motif'])
            feature['space'] |= set(codon['possible'])
        return feature

    def reference_for(self, id):
        if id not in self._reference_sequence:
            node = self.database['reference'].find_one({'head.id': id})
            if node:
                reference = Reference(self, node)
                self._reference_sequence[id] = reference
            else:
                self._reference_sequence[id] = None
        if id in self._reference_sequence:
            return self._reference_sequence[id]
        else:
            return None

    def accession_for(self, id):
        if id not in self._accession:
            node = self.database['accession'].find_one({'head.id': id})
            if node:
                self._accession[id] = node
            else:
                self._accession[id] = None
        if id in self._accession:
            return self._accession[id]
        else:
            return None

    def build_query(self, override, profile, kind='sample'):
        query = {}
        if profile is not None and profile in self.configuration['profile']:
            for k,v in self.configuration['profile'][profile][kind].items():
                query[k] = v
                
        if override:
            for k,v in override.items():
                query['head.{}'.format(k)] = v
        self.log.debug('Query is \n{}'.format(to_json(query)))
        return query

    def rebuild(self):
        existing_collections = self.database.collection_names()
        for table in self.configuration['table']:
            collection = self.database[table['collection']]
            if table['collection'] in existing_collections:
                existing_indexes = collection.index_information()
                for definition in table['index']:
                    if definition['name'] in existing_indexes:
                        self.log.info('dropping index %s on collection %s', definition['name'], table['collection'])
                        collection.drop_index(definition['name'])
                
            for definition in table['index']:
                self.log.info('building index %s on collection %s', definition['name'], table['collection'])
                collection.create_index(definition['key'], name=definition['name'], unique=definition['unique'])

    def describe(self, section):
        if section is not None:
            if section == '':
                pass
    

    def search_blat(self):
        result = None
        command = self.configuration['command']['blat']
        process = Popen(
            args=command['arguments'],
            cwd=command['cwd'],
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate(input=self.to_fasta().read().encode('utf8'))
        if output: result = StringIO(output.decode('utf8'))
        return result

    def populate_accession(self, format):
        def save_accession_record(collection, accession):
            if accession is not None:
                accession['organism name'] = 'mus musculus'
                if 'description' in accession:
                    accession['simple description'] = accession['description'].lower()
                    for pattern in self.configuration['strain']['mus musculus']:
                        if pattern['expression'].search(accession['simple description']):
                            accession['strain'] = pattern['name']
                            break
                        
                # fix the sequence structure
                if 'sequence' in accession:
                    if accession['sequence']:
                        accession['sequence'] = {
                            'nucleotide': accession['sequence'],
                            'strand': True,
                            'read frame': 0
                        }
                    else:
                        del accession['sequence']
                        
                if 'sequence' in accession:
                    record = { 'head': { 'id': accession['accession number'] }, 'body': accession }
                    for k in (
                        'strain',
                        'organism name',
                        'format'
                    ):
                        if k in accession:
                            record['head'][k] = accession[k]
                            del accession[k]
                    collection.save(record)
                else:
                    self.log.error('refusing to save accession record %s with missing sequence', accession['accession number'])

        if format in self.configuration['accession header']:
            parser = self.configuration['accession header'][format]
        else:
            self.log.critical('unknown fasta header format %s', format)
            raise SystemExit(1)
            
        collection = self.database['accession']
        accession = None
        for line in sys.stdin:
            if line is None: break
            else:
                line = line.strip()
                if line:
                    if line[0] == '>':
                        # at the start of a new record save the completed one
                        save_accession_record(collection, accession)
                        
                        # initialize a new accession
                        accession = { 'header': line, 'sequence': '', 'format': format }
                        match = parser.search(line)
                        if match:
                            parsed = dict(((k.replace('_', ' '),v) for k,v in match.groupdict().items()))
                            for k,v in parsed.items():
                                if v: accession[k] = v
                    else:
                        if self.configuration['expression']['nucleotide sequence'].search(line):
                            accession['sequence'] += line.upper()
        # save the last one
        save_accession_record(collection, accession)

    def populate_reference(self, path):
        def save_reference_record(collection, reference):
            if reference is not None:
                # fix the sequence structure
                if 'sequence' in reference:
                    if reference['sequence']:
                        reference['sequence'] = { 'nucleotide': reference['sequence'], 'strand': True }
                        
                        # fix the read frame
                        reference['sequence']['read frame'] = reference['read frame']
                        del reference['read frame']
                    else:
                        del reference['sequence']
                        
                if 'sequence' in reference:
                    record = { 'head': { 'id': reference['allele name'] }, 'body': reference }
                    for k in (
                        'framed',
                        'functionality',
                        'region',
                        'strain',
                        'organism name',
                        'confirmed genomoic dna',
                        'identified genomoic dna',
                        'accession'
                    ):
                        if k in reference:
                            record['head'][k] = reference[k]
                            del reference[k]
                    if record['head']['confirmed genomoic dna'] and record['head']['identified genomoic dna']:
                        record['head']['verified'] = True
                    else:
                        record['head']['verified'] = False
                    
                    if self.accession_for(record['head']['accession']) is None:
                        self.log.error('accession %s is missing', record['head']['accession'])
                    collection.save(record)
                else:
                    self.log.error('refusing to save reference record %s with missing sequence', reference['id'])

        collection = self.database['reference']
        reference = None
        with io.open(path, 'rb') as fasta:
            for line in fasta:
                if line:
                    line = line.strip().decode('utf8')
                    if line[0] == '>':
                        # at the start of a new record save the completed one
                        save_reference_record(collection, reference)
                            
                        # initialize a new reference
                        reference = {
                            'imgt meta': line,
                            'sequence': '',
                            'confirmed genomoic dna': True,
                            'identified genomoic dna': True,
                        }
                        match = self.configuration['expression']['imgt fasta header'].search(line)
                        if match:
                            parsed = dict(((k.replace('_', ' '),v) for k,v in match.groupdict().items()))
                            for k,v in parsed.items():
                                if v:
                                    try:
                                        if k in ['start', 'end', 'length']:
                                            v = int(v)
                                        elif k == 'read frame':
                                            if v == 'NR':
                                                v = None
                                            else:
                                                v = int(v)
                                        elif k == 'functionality':
                                            if '(' in v:
                                                reference['confirmed genomoic dna'] = False
                                                v = v.strip('()')
                                            elif '[' in v:
                                                reference['identified genomoic dna'] = False
                                                v = v.strip('[]')
                                            if v == 'ORF': v = 'O'
                                        elif k == 'organism name':
                                            v = v.lower()
                                            
                                    except ValueError as e:
                                        self.log.error('unable to decode %s as int for %s', v, k)
                                    else:
                                        if v: reference[k] = v
                                            
                            # fix the region for heavy chain
                            reference['region'] += 'H'
                            
                            # fix the accession strand
                            reference['accession strand'] = True
                            if 'polarity' in reference:
                                if reference['polarity'] == 'rev-compl':
                                    reference['accession strand'] = False
                                del reference['polarity']
                                    
                            # fix the start offset to zero based
                            if 'start' in reference: reference['start'] -= 1
                            
                            # fix the read frame offset to zero based
                            if 'read frame' in reference:
                                reference['read frame'] -= 1
                                reference['framed'] = True
                            else:
                                reference['framed'] = False
                                reference['read frame'] = 0
                    else:
                        if self.configuration['expression']['nucleotide sequence'].search(line):
                            reference['sequence'] += line.upper()
                else:
                    break
                    
            # save the last one if it has a sequence
            save_reference_record(collection, reference)

    def accession_to_fasta(self):
        buffer = StringIO()
        collection = self.database['accession']
        cursor = collection.find()
        for node in cursor:
            buffer.write(node['body']['header'])
            buffer.write('\n')
            sequence = Sequence(self,node['body']['sequence'])
            begin = 0
            end = 0
            while end < sequence.length:
                end = min(begin + self.configuration['constant']['fasta line length'], sequence.length)
                buffer.write(sequence.nucleotide[begin:end])
                buffer.write('\n')
                begin = end
        buffer.seek(0)
        print(buffer.read())

    def reference_to_blast_fasta(self, query, profile):
        q = self.build_query(query, profile, 'reference')
        buffer = StringIO()
        collection = self.database['reference']
        cursor = collection.find(q)
        for node in cursor:
            reference = Reference(self, node)
            buffer.write('>')
            buffer.write(reference.id)
            buffer.write('\n')
            begin = 0
            end = 0
            while end < reference.sequence.length:
                end = min(begin + self.configuration['constant']['fasta line length'], reference.sequence.length)
                buffer.write(reference.sequence.nucleotide[begin:end])
                buffer.write('\n')
                begin = end
        buffer.seek(0)
        print(buffer.read())

    def reference_to_auxiliary(self, query=None, profile='default'):
        q = self.build_query(query, profile, 'reference')
        buffer = StringIO()
        collection = self.database['reference']
        cursor = collection.find(q)
        for node in cursor:
            reference = Reference(self, node)
            if reference.framed:
                buffer.write(reference.id)
                buffer.write('\t')
                buffer.write(str(reference.sequence.read_frame))
                buffer.write('\t')
                buffer.write(reference.region)
                buffer.write('\n')
        buffer.seek(0)
        print(buffer.read())

    def drop_library(self, library):
        if library:
            collection = self.database['sample']
            try:
                result = collection.delete_many({ 'head.library': library })
                if result:
                    self.log.info('dropped %d samples from %s', result.deleted_count, library)
            except BulkWriteError as e:
                self.log.critical(e.details)

    def populate(self, library, strain, drop):
        count = 0
        if not library:
            raise ValueError('must specify a library to populate') 
        block = Block(self)
        collection = self.database['sample']
        if drop:
            self.drop_library(library)

        while block.fill(library, strain):
            try:
                result = collection.insert_many(block.document)
            except BulkWriteError as e:
                self.log.critical(e.details)
                raise SystemExit()
            count += block.size
            self.log.info('%s so far', count)

    def fasta(self, query=None, limit=None, skip=None, profile='default'):
        q = self.build_query(query, profile)
        collection = self.database['sample']
        cursor = collection.find(q)
        if limit is not None:
            cursor.limit(limit)
        if skip is not None:
            cursor.skip(skip)
            
        for node in cursor:
            sample = Sample(self, node)
            print(sample.fasta)
        cursor.close()

    def fastq(self, query=None, limit=None, skip=None, profile='default'):
        q = self.build_query(query, profile)
        collection = self.database['sample']
        cursor = collection.find(q)
        if limit is not None:
            cursor.limit(limit)
        if skip is not None:
            cursor.skip(skip)
            
        for node in cursor:
            sample = Sample(self, node)
            print(sample.fastq)
        cursor.close()

    def view(self, query=None, limit=None, skip=None, profile='default'):
        q = self.build_query(query, profile)
        collection = self.database['sample']
        cursor = collection.find(q)
        if limit is not None:
            cursor.limit(limit)
        if skip is not None:
            cursor.skip(skip)
            
        for node in cursor:
            sample = Sample(self, node)
            sample.view(profile)
        cursor.close()

    def info(self, query, limit=None, skip=None, profile='default'):
        q = self.build_query(query, profile)
        collection = self.database['sample']
        cursor = collection.find(q)
        if limit is not None:
            cursor.limit(limit)
        if skip is not None:
            cursor.skip(skip)
            
        for node in cursor:
            sample = Sample(self, node)
            sample.info(profile)
        cursor.close()

    def count(self, query=None, profile='default'):
        q = self.build_query(query, profile)
        collection = self.database['sample']
        print(collection.count(q))

    def simulate(self, library, strain, json, alignment, profile):
        block = Block(self)
        while block.fill(library, strain):
            block.simulate(json, alignment, profile)

    def execute(self, cmd):
        if cmd.action == 'accession-populate':
            self.populate_accession(cmd.instruction['format'])
            
        if cmd.action == 'accession-fasta':
            self.accession_to_fasta()
            
        if cmd.action == 'ref-populate':
            for path in cmd.instruction['path']:
                self.populate_reference(path)
                
        elif cmd.action == 'ref-fasta':
            self.reference_to_blast_fasta(cmd.query, cmd.instruction['profile'])
            
        elif cmd.action == 'ref-igblast-aux':
            self.reference_to_auxiliary(cmd.query, cmd.instruction['profile'])
                
        elif cmd.action == 'populate':
            self.populate(
                cmd.instruction['library'],
                cmd.instruction['strain'],
                cmd.instruction['drop'])
            
        elif cmd.action == 'drop':
            self.drop_library(cmd.instruction['library'])
            
        elif cmd.action == 'view':
            self.view(
                cmd.query,
                cmd.instruction['limit'],
                cmd.instruction['skip'],
                cmd.instruction['profile'])
            
        elif cmd.action == 'count':
            self.count(cmd.query, cmd.instruction['profile'])
            
        elif cmd.action == 'fasta':
            self.fasta(
                cmd.query,
                cmd.instruction['limit'],
                cmd.instruction['skip'],
                cmd.instruction['profile'])
            
        elif cmd.action == 'fastq':
            self.fastq(
                cmd.query,
                cmd.instruction['limit'],
                cmd.instruction['skip'],
                cmd.instruction['profile'])
            
        elif cmd.action == 'info':
            self.info(
                cmd.query,
                cmd.instruction['limit'],
                cmd.instruction['skip'],
                cmd.instruction['profile'])
            
        elif cmd.action == 'simulate':
            self.simulate(
                cmd.instruction['library'], 
                cmd.instruction['strain'], 
                cmd.instruction['json'], 
                cmd.instruction['alignment'],
                cmd.instruction['profile'])
            
        elif cmd.action == 'rebuild':
            self.rebuild()


def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    
    configuration['interface']['prototype']['profile']['parameter']['choices'] = list(configuration['profile'].keys())
    configuration['interface']['prototype']['profile']['parameter']['help'] = \
        '[ {} ]'.format(' | '.join(sorted(configuration['profile'].keys())))

    cmd = CommandLineParser(configuration['interface'])
    if cmd.sectioned and cmd.action is None:
        cmd.help()
    else:
        logging.getLogger().setLevel(log_levels[cmd.instruction['verbosity']])
        pipeline = Pipeline()
        try:
            pipeline.execute(cmd)

        except ValueError as e:
            logging.getLogger('main').critical(e)
            sys.exit(1)
        except(KeyboardInterrupt, SystemExit) as e:
            pipeline.close()
            sys.exit(1)
        
    pipeline.close()
    sys.exit(0)
            
            

if __name__ == '__main__':
    main()

