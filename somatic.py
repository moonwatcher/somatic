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
import io

from io import StringIO, BytesIO
from datetime import timedelta, datetime
from argparse import ArgumentParser
from subprocess import Popen, PIPE

from bson.objectid import ObjectId
from pymongo.son_manipulator import SONManipulator
from pymongo import MongoClient, DESCENDING, ASCENDING

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
    'yes/no question': ['Y', 'N'],
    'regions': ['VH', 'DH', 'JH'],
    'expression': {
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
            (?P<imgt_accession>[^|]*)\|
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
        'nucleotide sequence': re.compile('^[atcgnyATCGNY]+$')
    },
    'command': {
        'igblast command': [
            'igblastn',
            '-germline_db_V', 'database/mouse_imgt_vh',
            '-germline_db_J', 'database/mouse_imgt_jh',
            '-germline_db_D', 'database/mouse_imgt_dh',
            '-num_alignments_V', '5',
            '-num_alignments_J', '5',
            '-num_alignments_D', '5',
            '-organism', 'mouse',
            '-domain_system', 'imgt',
            '-query', '-',
            '-auxiliary_data', 'optional_file/mouse_gl.aux',
            '-show_translation',
            '-outfmt',
            '7 qseqid sseqid qstart qend sstart send gapopen gaps mismatch pident bitscore evalue length sstrand',
        ],
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
                    'metavar': 'NAME', 
                    'required': True
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
            'premature termination': {
                'flag': [
                    '-T', 
                    '--premature-termination'
                ], 
                'parameter': {
                    'choices': [
                        'Y', 
                        'N'
                    ], 
                    'dest': 'premature termination', 
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
                        'library'
                    ], 
                    'instruction': {
                        'description': 'match each read in file to regions with igblast and store results in the library. takes read from stdin', 
                        'help': 'populate samples for library', 
                        'name': 'populate'
                    }
                }, 
                {
                    'argument': [
                        'library', 
                        'id', 
                        'gapped', 
                        'valid', 
                        'in frame', 
                        'premature termination'
                    ], 
                    'instruction': {
                        'help': 'view sample alignment', 
                        'name': 'view'
                    }
                }, 
                {
                    'argument': [
                        'library', 
                        'id', 
                        'gapped', 
                        'valid', 
                        'in frame', 
                        'premature termination'
                    ], 
                    'instruction': {
                        'help': 'view JSON sample record', 
                        'name': 'info'
                    }
                }, 
                {
                    'argument': [
                        'library', 
                        'json', 
                        'alignment'
                    ], 
                    'instruction': {
                        'help': 'view JSON sample record', 
                        'name': 'simulate'
                    }
                }, 
                {
                    'argument': [
                        'path'
                    ], 
                    'instruction': {
                        'help': 'load imgt reference sequences from fasta file', 
                        'name': 'load-reference'
                    }
                }, 
                {
                    'argument': [
                        'region'
                    ], 
                    'instruction': {
                        'help': 'dump reference sequences to fasta', 
                        'name': 'to-blast-fasta'
                    }
                }, 
                {
                    'argument': [
                        'region'
                    ], 
                    'instruction': {
                        'help': 'dump reference sequences to igblast auxiliary file', 
                        'name': 'to-auxiliary'
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
    'complement': { 
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'R': 'Y',
        'Y': 'R',
        'K': 'M',
        'M': 'K',
        'S': 'S',
        'W': 'W',
        'B': 'V',
        'D': 'G',
        'H': 'C',
        'V': 'A',
        'N': 'N'
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
        'buffer size': 128,
        'fasta line length': 80,
        'minimum v identity': 0.7, 
        'minimum d identity': 0.7, 
        'minimum j identity': 0.8, 
        'minimum v alignment': 45, 
        'minimum d alignment': 4, 
        'minimum j alignment': 30, 
        'd overlap penalty factor': 0.5
    }
}

class CommandLineParser(object):
    def __init__(self, node):
        self.node = node
        self.parser = ArgumentParser(**self.node['instruction'])
        self.load()

    @property
    def sectioned(self):
        return 'section' in self.node and 'action' in self.node['section'] and self.node['section']['action']

    def load(self):
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

    def parse(self):
        instruction = {}
        arguments = vars(self.parser.parse_args())
        for k,v in arguments.items():
            if k is not None and v is not None:
                instruction[k] = v
                
        if self.sectioned and 'action' not in instruction:
            self.parser.print_help()
            return None
        else:
            return instruction

    def query(self, command):
        query = {}
        for k in [
            'library',
            'id'
        ]:
            if k in command:
                query[k] = command[k]
                
        for k in [
            'gapped',
            'valid',
            'in frame',
            'premature termination'
        ]:
            if k in command:
                if command[k] == 'Y':
                    query[k] = True
                elif command[k] == 'N':
                    query[k] = False
        return query


class Sequence(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Sequence')
        self.pipeline = pipeline 
        self.node = node
        self._reversed = None
        if self.node is None:
            self.node = {
                'read frame': 0,
                'strand': True
            }
        elif 'codon' in self.node:
            del self.node['codon']

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

    def reset(self):
        for k in ['codon']:
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
                    if acid in self.configuration['stop codon repertuar']:
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
            self.node = {}
        
        if 'sequence' in self.node and isinstance(self.node['sequence'], dict):
            self.node['sequence'] = Sequence(self.pipeline, self.node['sequence'])
        else:
            self.node['sequence'] = Sequence(self.pipeline)

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def id(self):
        return self.node['allele name']

    @property
    def sequence(self):
        return self.node['sequence']

    @property
    def region(self):
        return self.node['region']

    @property
    def strain(self):
        if 'strain' in self.node:
            return self.node['strain']
        else:
            return None

    @property
    def framed(self):
        return self.node['framed']

    @property
    def functionality(self):
        result = None
        if 'functionality' in self.node:
            if 'F' in self.node['functionality']:
                result = 'F'
            if 'P' in self.node['functionality']:
                 result = 'P'
            if 'ORF' in self.node['functionality']:
                result = 'O'
        return result


class Sample(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Sample')
        self.pipeline = pipeline
        self.node = node
        self.initialize()

    @property
    def configuration(self):
        return self.pipeline.configuration

    def initialize(self):
        def transform_hit(node):
            if 'query' in node and isinstance(node['query'], dict):
                node['query'] = Sequence(self.pipeline, node['query'])
                
            if 'subject' in node and isinstance(node['subject'], dict):
                node['subject'] = Sequence(self.pipeline, node['subject'])

        if self.node is None:
            self.node = {}
        
        if 'sequence' in self.node and isinstance(self.node['sequence'], dict):
            self.node['sequence'] = Sequence(self.pipeline, self.node['sequence'])
        else:
            self.node['sequence'] = Sequence(self.pipeline)
            
        if 'hit' in self.node:
            for hit in self.node['hit']:
                transform_hit(hit)
        else:
            self.node['hit'] = []
            
        if 'region' in self.node:
            for value in self.node['region'].values():
                transform_hit(value)
        else:
            self.node['region'] = {}

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
                return node.node
                
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
    def id(self):
        return self.node['id']

    @id.setter
    def id(self, value):
        self.node['id'] = value

    @property
    def strand(self):
        return self.node['strand']

    @strand.setter
    def strand(self, value):
        self.node['strand'] = value

    @property
    def library(self):
        return self.node['library']

    @library.setter
    def library(self, value):
        self.node['library'] = value

    @property
    def sequence(self):
        return self.node['sequence']

    @property
    def valid(self):
        return self.node['valid']

    @valid.setter
    def valid(self, value):
        self.node['valid'] = value

    @property
    def hit(self):
        return self.node['hit']

    @property
    def region(self):
        if 'region' not in self.node:
            self.node['region'] = {}
        return self.node['region']

    @property
    def framed(self):
        return self.node['framed']

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

    def invalidate(self, message):
        self.valid = False
        self.node['error'] = message
        self.log.warning('sample %s invalidated because %s', self.id, message)

    def invalidate_hit(self, hit, message):
        hit['valid'] = False
        hit['error'] = message
        self.log.warning('%dbp %s hit to %s on %s invalidated because %s', 
        hit['alignment length'], hit['region'], hit['subject id'], self.id, message)

    def _check_hit_strand(self, reference, hit):
        # only DH regions are allowed to align to the oposite strand
        if hit['region'] != 'DH' and hit['subject strand'] != reference.sequence.strand:
            self.invalidate_hit(hit, 'hit aligns to the wrong strand')
        return hit['valid']

    def _set_frame_and_strand(self, hit):
        self.node['framed'] = True
        self.sequence.read_frame = (hit['query start'] + hit['subject'].read_frame) % 3
        self.strand = hit['subject strand']
        self._assign_frame()
        self.log.debug('frame set for %s by %s', self.id, hit['subject id'])

    def _default_analysis_values(self):
        for hit in self.hit: hit['valid'] = True
        self.node['framed'] = False
        self.node['valid'] = True
        self.node['premature termination'] = False
        self.node['in frame'] = False

    def _check_for_gapped_alignment(self):
        for hit in self.hit:
            if hit['gap openings'] > 0:
                hit['gapped'] = True
                self.invalidate_hit(hit, 'hit contains a gapped alignment')

    def _complement_from_reference(self):
        for hit in self.hit:
            hit['picked'] = False
            hit['query start'] -= 1
            hit['query'] = self.sequence.crop(hit['query start'], hit['query end'])
            
            reference = self.pipeline.reference_for(hit['subject id'])
            
            if hit['subject strand'] == reference.sequence.strand:
                hit['subject start'] -= 1
                hit['subject'] = reference.sequence.crop(hit['subject start'], hit['subject end'])
            else:
                hit['subject start'] = reference.sequence.length - hit['subject start']
                hit['subject end'] = reference.sequence.length - hit['subject end'] + 1
                hit['subject'] = reference.sequence.reversed.crop(hit['subject start'], hit['subject end'])

            hit['framed'] = reference.framed
            hit['functionality'] = reference.functionality
            if reference.strain:
                hit['strain'] = reference.strain
                
            self._check_hit_strand(reference, hit)

    def _pick_jh_region(self):
        picked = False
        for hit in self.hit:
            if (hit['valid'] and
                hit['framed'] and 
                hit['region'] == 'JH' and 
                hit['identical'] >=  self.configuration['constant']['minimum j identity'] and
                hit['alignment length'] >= self.configuration['constant']['minimum j alignment']):
                
                self.region['JH'] = hit
                hit['picked'] = True
                
                # JH match sets the read frame and strand for the sample
                if not self.node['framed']: self._set_frame_and_strand(hit)
                picked = True
                break
        # if not picked: self.invalidate('no satisfactory JH region match found')
        return picked

    def _assign_frame(self):
        for hit in self.hit:
            if hit['valid']:
                hit['query'].read_frame = 2 - (hit['query start'] - self.sequence.read_frame - 1) % 3
                if hit['framed']:
                    if hit['query'].read_frame == hit['subject'].read_frame:
                        hit['in frame'] = 'Y'
                    else:
                        hit['in frame'] = 'N'
                else:
                    hit['in frame'] = '*'

    def _pick_vh_region(self):
        picked = False
        for hit in self.hit:
            if (hit['valid'] and
                hit['framed'] and 
                hit['region'] == 'VH' and 
                hit['identical'] >=  self.configuration['constant']['minimum v identity'] and
                hit['alignment length'] >= self.configuration['constant']['minimum v alignment']):
                
                self.region['VH'] = hit
                hit['picked'] = True

                # if JH didn't set the frame and strand a VH match will
                if not self.node['framed']: self._set_frame_and_strand(hit)
                picked = True
                break
        # if not picked: self.invalidate('no satisfactory VH region match found')
        return picked

    def _pick_dh_region(self):
        picked = False
        possible = []
        for hit in self.hit:
            if hit['valid'] and hit['region'] == 'DH':
                h = {}
                for k,v in hit.items():
                    if k not in ('query', 'subject'): h[k] = v
                reference = self.pipeline.reference_for(hit['subject id'])
                    
                h['overlap'] = 0
                # check for overlap with VH
                if 'VH' in self.region and self.region['VH']['query end'] > hit['query start']:
                    h['query start'] = self.region['VH']['query end']
                    overlap = h['query start'] - hit['query start']
                    h['subject start'] += overlap
                    h['overlap'] += overlap
                    
                # check for overlap with JH
                if 'JH' in self.region and hit['query end'] > self.region['JH']['query start']:
                    h['query end'] = self.region['JH']['query start']
                    overlap = hit['query end'] - h['query end']
                    h['subject end'] -= overlap
                    h['overlap'] += overlap
                    
                # correct the score for the overlap and trim the sequences
                h['score'] = hit['bit score'] - h['overlap'] * self.configuration['d overlap penalty factor']
                h['alignment length'] = h['query end'] - h['query start']
                h['query'] = self.sequence.crop(h['query start'], h['query end'])

                if h['subject strand'] == reference.sequence.strand:
                    h['subject'] = reference.sequence.crop(h['subject start'], h['subject end'])
                else:
                    h['subject'] = reference.sequence.reversed.crop(h['subject start'], h['subject end'])
                                    
                # filter DH hits that match the criteria
                if h['query end'] > h['query start']:
                    similar = 0
                    different = 0
                    longest = 0
                    stretch = 0
                    for index in range(h['query'].length):
                        if h['subject'].nucleotide[index] == h['query'].nucleotide[index]:
                            similar += 1
                            stretch += 1
                            longest = max(longest, stretch)
                        else:
                            different += 1
                            stretch = 0
                    h['identical'] = float(similar) / float(h['alignment length'])
                    if (longest >= self.configuration['constant']['minimum d alignment'] and
                        h['identical'] >= self.configuration['constant']['minimum d identity']):
                        possible.append(h)
                        
        # pick the highest preforming match
        if possible:
            possible.sort(reverse=True, key=lambda x: x['score'])
            self.region['DH'] = possible[0]
            for h in self.hit:
                if h['region'] == 'DH' and h['subject id'] == self.region['DH']['subject id']:
                    h['picked'] = True
            if not self.node['framed']: self._set_frame_and_strand(self.region['DH'])
            picked = True
            
        if not picked: self.log.debug('no DH region match passed quality assessment for %s', self.id)
            
        return picked

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
            for index,codon in enumerate(self.region['JH']['query'].codon):
                if codon == 'W':
                    offset = self.region['JH']['query'].read_frame + (index + 1) * 3
                    break
                    
            if offset is not None:
                end = self.region['JH']['query start'] + offset
                
            if start and end:
                self.region['CDR3'] = {
                    'query start': start,
                    'query end': end,
                    'region': 'CDR3',
                    'query': self.sequence.crop(start, end)
                }

    def _check_for_stop_codon(self):
        if '*' in self.sequence.codon:
            self.node['premature termination'] = True
        else:
            self.node['premature termination'] = False
        return not self.node['premature termination']

    def _check_for_aligned_frames(self):
        if 'JH' in self.region and 'VH' in self.region:
            if self.region['VH']['in frame'] == 'Y' and self.region['JH']['in frame'] == 'Y':
                self.node['in frame'] = True
            else:
                self.node['in frame'] = False
        else:
            self.node['in frame'] = False

    def _identify_v_d_junction(self):
        if 'DH' in self.region and 'VH' in self.region:
            if self.region['VH']['query end'] < self.region['DH']['query start']:
                start = self.region['VH']['query end']
                end = self.region['DH']['query start']
                self.region['VD junction'] = {
                    'query start': start,
                    'query end': end,
                    'region': 'VD junction',
                    'query': self.sequence.crop(start, end)
                }
        return True

    def _identify_d_j_junction(self):
        if 'DH' in self.region and 'JH' in self.region:
            if self.region['JH']['query start'] > self.region['DH']['query end']:
                start = self.region['DH']['query end']
                end = self.region['JH']['query start']
                self.region['DJ junction'] = {
                    'query start': start,
                    'query end': end,
                    'region': 'VD junction',
                    'query': self.sequence.crop(start, end)
                }
        return True

    def _identify_v_j_junction(self):
        if 'DH' not in self.region and 'JH' in self.region and 'VH' in self.region:
            if self.region['JH']['query start'] > self.region['VH']['query end']:
                start = self.region['VH']['query end']
                end = self.region['JH']['query start']
                self.region['VJ junction'] = {
                    'query start': start,
                    'query end': end,
                    'region': 'VJ junction',
                    'query': self.sequence.crop(start, end)
                }
        return True

    def analyze(self):
        self._default_analysis_values()
        if self.valid: self._check_for_gapped_alignment()
        if self.valid: self._complement_from_reference()
        if self.valid: self._pick_jh_region()
        if self.valid: self._pick_vh_region()
        if self.valid: self._check_for_aligned_frames()
        if self.valid: self._pick_dh_region()
        if self.valid: self._identify_v_j_junction()
        if self.valid: self._identify_v_d_junction()
        if self.valid: self._identify_d_j_junction()
        if self.valid: self._check_for_stop_codon()
        if self.valid: self._identify_cdr3()

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
            
        return json.dumps(self.node, sort_keys=False, ensure_ascii=False, indent=4, default=handler)

    def to_summary(self):
        if self.hit:
            buffer = StringIO()
            buffer.write(self.id)
            buffer.write(' | ')
            for r in ('VH', 'DH', 'JH'):
                if r in self.region:
                    buffer.write(r)
                    buffer.write(' : ')
                    buffer.write(self.region[r]['subject id'])
                    buffer.write(' | ')
            
            if self.node['in frame']:
                buffer.write('in frame')
            else:
                buffer.write('out of frame')
                
            if self.node['premature termination']:
                buffer.write(' | prematurely terminated')
                
            buffer.seek(0)
            return buffer.read()
        else:
            return ''

    def to_alignment_diagram(self):
        def print_gap(buffer, framed):
            if framed: buffer.write(' ')

        if self.hit:
            buffer = StringIO()
            offset = {
                'gene': len('gene'),
                'strain': len('strain')
            }
            
            c = [ len(hit['subject id']) for hit in self.hit if hit['valid'] ]
            if c: offset['gene'] = max(max(c), offset['gene'])
            
            c = [ len(s) for s in self.pipeline.strains ]
            offset['strain'] = max(max(c), offset['strain'])
            
            offset['alignment'] = offset['gene'] + offset['strain'] + 15
            offset['sample frame'] = self.sequence.read_frame
            
            # print a summary
            buffer.write('{: <{}}{}'.format('summary', offset['alignment'], self.to_summary()))
            buffer.write('\n')
            
            # print a title
            buffer.write('{: <{}} {: <{}} {: <4} {} {} {} {} '.format('gene', offset['gene'], 'strain', offset['strain'], 'R', 'F', 'I', 'P', 'S'))
            
            # print the sample coordinate system
            if offset['sample frame'] > 0:
                buffer.write('{: <{}}'.format(0, offset['sample frame']))
                print_gap(buffer, self.framed)
            
            gap = 3 if self.framed else 6
            for index in range(offset['sample frame'], self.sequence.length, gap):
                buffer.write('{: <{}}'.format(index, gap))
                print_gap(buffer, self.framed)
            buffer.write('\n')
            
            # print the sample nucleotide sequence
            buffer.write(' ' * offset['alignment'])
            if offset['sample frame'] > 0:
                buffer.write(self.sequence.nucleotide[0:offset['sample frame']])
                print_gap(buffer, self.framed)
            for index in range(offset['sample frame'], self.sequence.length, 3):
                buffer.write(self.sequence.nucleotide[index:index + 3])
                print_gap(buffer, self.framed)
            buffer.write('\n')
            
            # print the sample quality sequence
            buffer.write(' ' * offset['alignment'])
            if offset['sample frame'] > 0:
                buffer.write(self.sequence.quality[0:offset['sample frame']])
                print_gap(buffer, self.framed)
            for index in range(offset['sample frame'], self.sequence.length, 3):
                buffer.write(self.sequence.quality[index:index + 3])
                print_gap(buffer, self.framed)
            buffer.write('\n')
            
            if self.framed:
                # print the sample codon sequence
                buffer.write(' ' * offset['alignment'])
                if offset['sample frame'] > 0:
                    buffer.write(' ' * offset['sample frame'])
                    print_gap(buffer, self.framed)
                for codon in self.sequence.codon:
                    buffer.write('{: <3}'.format(codon))
                    print_gap(buffer, self.framed)
                buffer.write('\n')
            
            for hit in self.hit:
                if hit['valid']:
                    hit_offset = 0
                    if hit['query start'] > 0:
                        hit_offset = hit['query start']
                        if self.framed:
                            if offset['sample frame'] > 0:
                                hit_offset += int((hit['query start'] - offset['sample frame']) / 3) + 1
                            else:
                                hit_offset += (int(hit['query start'] / 3))
                            
                    buffer.write('{: <{}} {: <{}} {: <4} {} {} {} {} '.format(
                        hit['subject id'], offset['gene'], 
                        '' if 'strain' not in hit else hit['strain'], offset['strain'], 
                        hit['region'], 
                        hit['functionality'],
                        hit['in frame'] if 'in frame' in hit else ' ',
                        'Y' if hit['picked'] else 'N',
                        '+' if hit['subject strand'] else '-'))
                    
                    subject = hit['subject'].clone()
                    subject.read_frame = hit['query'].read_frame
                    
                    # print the hit nucleotide mask
                    mask = []
                    for index,n in enumerate(subject.nucleotide):
                        mask.append('-' if n == hit['query'].nucleotide[index] else n)
                    mask = ''.join(mask)
                        
                    if hit['query start'] > 0: buffer.write(' ' * hit_offset)
                    if subject.read_frame > 0:
                        buffer.write(mask[0:subject.read_frame])
                        print_gap(buffer, self.framed)
                        
                    for index in range(subject.read_frame, subject.length, 3):
                        buffer.write(mask[index:index + 3])
                        print_gap(buffer, self.framed)
                    buffer.write('\n')
                    
                    if self.framed:
                        # print the codon mask
                        display = False
                        mask = []
                        for index,c in enumerate(subject.codon):
                            if c == hit['query'].codon[index]: mask.append('-')
                            else:
                                mask.append(c)
                                display = True
                        if display:
                            buffer.write(' ' * offset['alignment'])
                            if hit['query start'] > 0: buffer.write(' ' * hit_offset)
                                
                            if subject.read_frame > 0:
                                buffer.write(' ' * subject.read_frame)
                                print_gap(buffer, self.framed)
                                
                            for codon in mask:
                                buffer.write('{: <3}'.format(codon))
                                print_gap(buffer, self.framed)
                            buffer.write('\n')
            buffer.seek(0)
            return buffer.read()
        else:
            return ''

    def view(self):
        print(self.to_alignment_diagram())

    def info(self):
        print(self.to_json())

    def __str__(self):
        return self.to_alignment_diagram()


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
        self.buffer.append(sample)
        self.lookup[sample.id] = sample

    def fill(self, library):
        if self.read():
            buffer = self.search()
            if buffer:
                self.parse_igblast(buffer)
                for sample in self.buffer:
                    sample.library = library
                    sample.analyze()
        return not self.empty

    def read(self):
        self.reset()
        state = None
        sample = Sample(self.pipeline)
        for line in sys.stdin:
            if line:
                line = line.strip()
                if line:
                    if state == None:
                        if line[0] == '@':
                            sample.id = line
                            state = 'id'
                            
                    elif state == 'id':
                        sample.sequence.nucleotide = line
                        state = 'sequence'
                        
                    elif state == 'sequence':
                        if line[0] == '+':
                            state = 'redundant'
                            
                    elif state == 'redundant':
                        sample.sequence.quality = line
                        self.add(sample)
                        state = None
                        if self.full:
                            break
                        else:
                            sample = Sample(self.pipeline)
            else:
                break
                
        return not self.empty

    def search(self):
        result = None
        process = Popen(
            args=self.configuration['command']['igblast command'],
            cwd='/Users/lg/code/somatic/igblast',
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate(input=self.to_fasta().read().encode('utf8'))
        if output:
            result = StringIO(output.decode('utf8'))
        return result

    def to_fasta(self):
        buffer = StringIO()
        for sample in self.buffer:
            buffer.write('>{}\n'.format(sample.id))
            begin = 0
            end = 0
            while end < sample.sequence.length:
                end = min(begin + self.configuration['constant']['fasta line length'], sample.sequence.length)
                buffer.write(sample.sequence.nucleotide[begin:end])
                buffer.write('\n')
                begin = end
        buffer.seek(0)
        return buffer

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
                    sample.hit.append(hit)

    def view(self):
        for sample in self.buffer:
            sample.view()

    def info(self):
        for sample in self.buffer:
            sample.info()

    def simulate(self, json, alignment):
        for sample in self.buffer:
            if json: sample.info()
            if alignment: sample.view()


class Pipeline(object):
    def __init__(self):
        self.log = logging.getLogger('Pipeline')
        self.configuration = configuration
        self.count = 0
        self._connection = None
        self._reference_sequence = {}
        self._stop_codon_feature = None
        self._strains = None
        
        self.configuration['complement'] = {}
        for k,v in self.configuration['iupac nucleic acid notation'].items():
            self.configuration['complement'][k] = v['reverse']
        
        stop = [ k for k,v in self.configuration['nucleic to amino'].items() if v == '*' ]
        motif = self.motif_for(stop)
        self.configuration['stop codon repertuar'] = motif['space']

    @property
    def strains(self):
        if self._strains == None:
            self._strains = list(self.database['reference'].distinct('strain'))
        return self._strains

    @property
    def connection(self):
        if self._connection is None:
            try:
                self._connection = MongoClient('mongodb://localhost')
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
            node = self.database['reference'].find_one({'allele name': id})
            if node:
                reference = Reference(self, node)
                self._reference_sequence[id] = reference
            else:
                self._reference_sequence[id] = None
        return self._reference_sequence[id]

    def rebuild(self):
        index = [
            {
                'collection': 'reference',
                'index': [
                    { 'key': [( 'allele name', ASCENDING )], 'unique': True, 'name': 'reference allele name' },
                    { 'key': [( 'region', ASCENDING )], 'unique': False, 'name': 'reference region' },
                    { 'key': [( 'subgroup', ASCENDING )], 'unique': False, 'name': 'reference subgroup' },
                    { 'key': [( 'gene name', ASCENDING )], 'unique': False, 'name': 'reference gene name' }
                ]
            },
            {
                'collection': 'sample',
                'index': [
                    { 'key': [( 'id', ASCENDING )], 'unique': True, 'name': 'sample id' },
                    { 'key': [( 'in frame', ASCENDING )], 'unique': False, 'name': 'sample in frame' },
                    { 'key': [( 'premature termination', ASCENDING )], 'unique': False, 'name': 'sample premature termination' },
                    { 'key': [( 'library', ASCENDING )], 'unique': False, 'name': 'sample library' },
                    { 'key': [( 'gapped', ASCENDING )], 'unique': False, 'name': 'sample gapped' },
                    { 'key': [( 'valid', ASCENDING )], 'unique': False, 'name': 'sample valid' },
                ]
            }
        ]
        existing_collections = self.database.collection_names()
        for table in index:
            collection = self.database[table['collection']]
            if table['collection'] in existing_collections:
                existing_indexes = collection.index_information()
                for definition in table['index']:
                    if definition['name'] in existing_indexes:
                        self.log.info('dropping index %s on collection %s', definition['name'], table['collection'])
                        collection.drop_index(definition['name'])
                
            for definition in table['index']:
                self.log.info('rebuilding index %s on collection %s', definition['name'], table['collection'])
                collection.create_index(definition['key'], name=definition['name'], unique=definition['unique'])

    def populate_reference(self, path):
        def save_reference_sample(collection, sample):
            if sample is not None:
                
                # fix the sequence structure
                if 'sequence' in sample:
                    if sample['sequence']:
                        sample['sequence'] = {
                            'nucleotide': sample['sequence'],
                            'strand': True
                        }
                        
                        # fix the read frame
                        if 'read frame' in sample:
                            sample['sequence']['read frame'] = sample['read frame']
                            del sample['read frame']
                            sample['framed'] = True
                        else:
                            sample['sequence']['read frame'] = 0
                            sample['framed'] = False
                    else:
                        del sample['sequence']
                        
                if 'sequence' in sample:
                    collection.save(sample)
                else:
                    self.log.error('refusing to save reference record %s with missing sequence', sample['id'])

        collection = self.database['reference']
        sample = None
        with io.open(path, 'rb') as fasta:
            for line in fasta:
                if line:
                    line = line.strip().decode('utf8')
                    if line[0] == '>':
                        # at the start of a new record save the completed one
                        save_reference_sample(collection, sample)
                            
                        # initialize a new sample
                        sample = { 'id': line, 'sequence': '' }
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
                                    except ValueError as e:
                                        self.log.error('unable to decode %s as int for %s', v, k)
                                    else:
                                        if v is not None:
                                            sample[k] = v
                                            
                            # fix the region for heavy chain
                            sample['region'] += 'H'
                            
                            # fix the strand
                            sample['strand'] = True
                            if 'polarity' in sample:
                                if sample['polarity'] == 'rev-compl': sample['strand'] = False
                                del sample['polarity']
                                    
                            # fix the start offset
                            if 'start' in sample: sample['start'] -= 1
                            
                            # fix the read frame offset to zero based
                            if 'read frame' in sample: sample['read frame'] -= 1
                    else:
                        if self.configuration['expression']['nucleotide sequence'].search(line):
                            sample['sequence'] += line.upper()
                else:
                    break
                    
            # save the last one if it has a sequence
            save_reference_sample(collection, sample)

    def reference_to_blast_fasta(self, region):
        buffer = StringIO()
        collection = self.database['reference']
        query = { 'region': region }
        cursor = collection.find(query)
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

    def reference_to_auxiliary(self, region):
        buffer = StringIO()
        collection = self.database['reference']
        query = { 'region': region }
        cursor = collection.find(query)
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

    def populate(self, library):
        block = Block(self)
        collection = self.database['sample']
        while block.fill(library):
            result = collection.insert_many(block.document)
            self.count += block.size
            self.log.debug('%s so far', self.count)

    def view(self, query):
        collection = self.database['sample']
        cursor = collection.find(query)
        for node in cursor:
            sample = Sample(self, node)
            sample.view()

    def info(self, query):
        collection = self.database['sample']
        cursor = collection.find(query)
        for node in cursor:
            sample = Sample(self, node)
            sample.info()

    def simulate(self, library, json, alignment):
        block = Block(self)
        while block.fill(library):
            block.simulate(json, alignment)


def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    
    cli = CommandLineParser(configuration['interface'])
    command = cli.parse()
    action = command['action']
    logging.getLogger().setLevel(log_levels[command['verbosity']])
    
    pipeline = Pipeline()
    if action == 'load-reference':
        for path in command['path']:
            pipeline.populate_reference(path)
            
    elif action == 'to-blast-fasta':
        pipeline.reference_to_blast_fasta(command['region'])
        
    elif action == 'to-auxiliary':
        pipeline.reference_to_auxiliary(command['region'])
            
    elif action == 'populate':
        pipeline.populate(command['library'])
        
    elif action == 'view':
        pipeline.view(cli.query(command))
        
    elif action == 'info':
        pipeline.info(cli.query(command))
        
    elif action == 'simulate':
        pipeline.simulate(command['library'], command['json'], command['alignment'])
        
    elif action == 'rebuild':
        pipeline.rebuild()
        
    pipeline.close()

if __name__ == '__main__':
    main()

