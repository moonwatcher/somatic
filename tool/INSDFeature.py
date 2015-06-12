#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import logging
import re
import json
import uuid
import io
import hashlib
import math

from io import StringIO, BytesIO
from datetime import timedelta, datetime
from argparse import ArgumentParser
from subprocess import Popen, PIPE

from bson.objectid import ObjectId
from pymongo.son_manipulator import SONManipulator
from pymongo import MongoClient, DESCENDING, ASCENDING
from pymongo.errors import BulkWriteError

import xmltodict
import urllib.request, urllib.parse, urllib.error
import urllib.parse
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from http.client import BadStatusLine

def to_json(node):
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4)

INSDFeature = [
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1",
                    "INSDInterval_to": "2500000"
                }
            }
        ],
        "INSDFeature_key": "source",
        "INSDFeature_location": "1..2500000",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "organism",
                        "INSDQualifier_value": "Mus musculus"
                    },
                    {
                        "INSDQualifier_name": "mol_type",
                        "INSDQualifier_value": "genomic DNA"
                    },
                    {
                        "INSDQualifier_name": "strain",
                        "INSDQualifier_value": "C57BL/6"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "taxon:10090"
                    },
                    {
                        "INSDQualifier_name": "chromosome",
                        "INSDQualifier_value": "12"
                    },
                    {
                        "INSDQualifier_name": "germline"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1",
                    "INSDInterval_to": "2500000"
                }
            }
        ],
        "INSDFeature_key": "misc_feature",
        "INSDFeature_location": "1..2500000",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "immunoglobulin heavy chain variable region, complete region containing V segments"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "300",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "201"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "complement(<201..300)",
        "INSDFeature_partial3": {
            "@value": "true"
        },
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Vipr2"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "300",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "201"
                }
            }
        ],
        "INSDFeature_key": "exon",
        "INSDFeature_location": "complement(201..300)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Vipr2"
                    },
                    {
                        "INSDQualifier_name": "number",
                        "INSDQualifier_value": "1"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "248",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "201"
                }
            }
        ],
        "INSDFeature_key": "CDS",
        "INSDFeature_location": "complement(<201..248)",
        "INSDFeature_partial3": {
            "@value": "true"
        },
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Vipr2"
                    },
                    {
                        "INSDQualifier_name": "inference",
                        "INSDQualifier_value": "similar to RNA sequence, mRNA:U09631"
                    },
                    {
                        "INSDQualifier_name": "codon_start",
                        "INSDQualifier_value": "1"
                    },
                    {
                        "INSDQualifier_name": "transl_table",
                        "INSDQualifier_value": "1"
                    },
                    {
                        "INSDQualifier_name": "product",
                        "INSDQualifier_value": "vasoactive intestinal peptide receptor 2"
                    },
                    {
                        "INSDQualifier_name": "protein_id",
                        "INSDQualifier_value": "CAJ55824.1"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "GI:90704842"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "GOA:Q1WWJ6"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "UniProtKB/TrEMBL:Q1WWJ6"
                    },
                    {
                        "INSDQualifier_name": "translation",
                        "INSDQualifier_value": "MRASVVLTCYCWLLVR"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "24522",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "11702"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "complement(11702..24522)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Zfp386"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "24522",
                        "INSDInterval_iscomp": {
                            "@value": "true"
                        },
                        "INSDInterval_to": "24433"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "17577",
                        "INSDInterval_iscomp": {
                            "@value": "true"
                        },
                        "INSDInterval_to": "17346"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "13394",
                        "INSDInterval_iscomp": {
                            "@value": "true"
                        },
                        "INSDInterval_to": "11702"
                    }
                ]
            }
        ],
        "INSDFeature_key": "mRNA",
        "INSDFeature_location": "complement(join(11702..13394,17346..17577,24433..24522))",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Zfp386"
                    },
                    {
                        "INSDQualifier_name": "product",
                        "INSDQualifier_value": "zinc finger protein 386 (Kruppel-like) isoform a"
                    },
                    {
                        "INSDQualifier_name": "transcription",
                        "INSDQualifier_value": "GAAAATGTGAGCGTGTTTTCCAGCTTTTGGAAGCCCCGGTGGCCCTTACGAAATCACGGGTGCTGGAGATGCAGGAGAGCCATGGAGATGGAGCAGTTGACCTTCCGCGACGTGGCCGTTGACTTCTCTCCGGATGAGTGGGAATGCCTGGACCCTCCTCAGCAGAGGCTGTATAGGGATGTGATGGTAGAGAACTACAGAAACCTGGTCTCTGTGGGTGAGGACAGCATCCCTGCAGAGCTCCCAACGCACTGCCACCATTTCATGCCTTTCCTTTCAAAAATGTCTGCCATAGCCTGTGCTTTTCATAAGCAGTACAAACCTCTTTCTCATCATTGTACCCAGAACTTTTCACCAGAGCAGATCCTGAAATATTCCTTCACAGAAGTGACTAATGGAAGTTATAGCCTTCAGCATTTATACTTAAGGGAAGATTGTGCAAATGTAGGTGGGAATAAAGCAGCAGCCTGCCATGGTAGACATAAGCAATGTCTGATCACAAGCAATACCAAAGACAACAGAAATCAAGAACAGAAAACAGCTCTGAAAAAACCTCATGCAAGGACAGATACTTTTGATGAGCCATGTGTCTATAACTGTAAATATCCACAGGAAATTTTTCAACATGACTTAACTTTGAAAGAAAGTTGGGAAGATCTGATAGATAGATTAATTTTTCACTCAAGCAATTTCATAGACAAAGAATTAAGAAACAGAGGACCGATCTCTAAAAGTGATCAAGTCGAGTCTTCCTATCCTGGAAGAGCCATACTCTTCAGTCAGCAAATACCAGTCTCCTGTGGCCAAAAGCACAATACTAATGATTGTGGATTGCCTGCTACTGATCCATTATTGCATAGCCAAAATCTAGATACAGATATTTGGGAAACCCCTTACATGTATAAGGAAATTAGCAAAGGTTTTAGTAGTGGATCTTTACTAAATAGTTGCAATGATGTGGTTATCTTAGAGAAATCTGGTCAGCGTGACAAAGTACAGAATGACTTTGACCATTGCTTGAGTCCCAACAACCATCAATGTAATCTTCCAAAGAAGCTTTTCAATTGTGATAATTTCTTTACTCAATGTTCAAAGCTTACTATCCAGCAGTGTAACTATATTCAAGATAATCATTACAAGAGTATAGATTGTGACACAATGTTTAATAAAACTTTAAACGTTACTAGGCGTAAAATTCAGCATTCCAAGAAATCCTACAAATGTAATGAATGTGGCAAAGCCTTTAAATATTGCTCAAGTTACAGAAAACACAGCATAATTCACACAGGAGAGAAACCCTACAAATGTAAATTATGTGGAAAATCCTTTACCCAGTGTGCAAGCCTTAAAAAACATCAGAGAATCCACACTGGAGAGAAACCCTACAGGTGTGAAGAATGTGGCAGAAGCTTTAACCACTACTCAATTCTTGGTCAGCACCAGAGGATTCACACAGGAGAGAAACCCTACAAGTGTAAACAGTGTGGCAAGTCCTTTACTCAGTGTTCAAGCCTTCAGAAACACCAAGTAATCCACACTGGAGAGAAGCCCTACAGATGTGCAGAATGTGGCAAATCTTTTACCCAAAACTCAACACTTAGTCAACACCAGAGAATTCATACTGGAGAGAAACCTTACAAATGTGAAGAGTGTGGAAAGGCCTTTACCCAATGTTCAAGCCTTAGAAAACATCAGAGAATTCACACTGGAGAGAAACCCTACAAATGTGAAGTATGTGGCAGAGCCTTTAACTGCCGCTCATCTTTTACTAAACACAAGAGAATTCACACAGGAGAGAAACCTTACAAATGTAAAGACTGTGACAAAGCTTTTATCCACTGCACAAATCTTATTCAACACCAGAGAATCCATACAGGAGAGAAACCATACAAATGTAATGAATGTGGGAAATCCTTTAGTCAATGTTCAAACCTTAGAAAGCATGAGAGAATTCATACATGAAAAAAACTAGCAAATGTCAAGAATCTAGTAAATCTTTGATCCATAGTCAAAAA"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "13394",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "11702"
                }
            }
        ],
        "INSDFeature_key": "exon",
        "INSDFeature_location": "complement(11702..13394)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Zfp386"
                    },
                    {
                        "INSDQualifier_name": "number",
                        "INSDQualifier_value": "3"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "24441",
                        "INSDInterval_iscomp": {
                            "@value": "true"
                        },
                        "INSDInterval_to": "24433"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "17577",
                        "INSDInterval_iscomp": {
                            "@value": "true"
                        },
                        "INSDInterval_to": "17346"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "13394",
                        "INSDInterval_iscomp": {
                            "@value": "true"
                        },
                        "INSDInterval_to": "11755"
                    }
                ]
            }
        ],
        "INSDFeature_key": "CDS",
        "INSDFeature_location": "complement(join(11755..13394,17346..17577,24433..24441))",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Zfp386"
                    },
                    {
                        "INSDQualifier_name": "inference",
                        "INSDQualifier_value": "similar to RNA sequence, mRNA:BC004747"
                    },
                    {
                        "INSDQualifier_name": "codon_start",
                        "INSDQualifier_value": "1"
                    },
                    {
                        "INSDQualifier_name": "transl_table",
                        "INSDQualifier_value": "1"
                    },
                    {
                        "INSDQualifier_name": "product",
                        "INSDQualifier_value": "zinc finger protein 386 (Kruppel-like) isoform a"
                    },
                    {
                        "INSDQualifier_name": "protein_id",
                        "INSDQualifier_value": "CAJ55825.1"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "GI:90704843"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "GOA:Q1WWJ5"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "InterPro:IPR001909"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "InterPro:IPR007087"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "InterPro:IPR013087"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "InterPro:IPR015880"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "MGI:MGI:1930708"
                    },
                    {
                        "INSDQualifier_name": "db_xref",
                        "INSDQualifier_value": "UniProtKB/TrEMBL:Q1WWJ5"
                    },
                    {
                        "INSDQualifier_name": "translation",
                        "INSDQualifier_value": "MEMEQLTFRDVAVDFSPDEWECLDPPQQRLYRDVMVENYRNLVSVGEDSIPAELPTHCHHFMPFLSKMSAIACAFHKQYKPLSHHCTQNFSPEQILKYSFTEVTNGSYSLQHLYLREDCANVGGNKAAACHGRHKQCLITSNTKDNRNQEQKTALKKPHARTDTFDEPCVYNCKYPQEIFQHDLTLKESWEDLIDRLIFHSSNFIDKELRNRGPISKSDQVESSYPGRAILFSQQIPVSCGQKHNTNDCGLPATDPLLHSQNLDTDIWETPYMYKEISKGFSSGSLLNSCNDVVILEKSGQRDKVQNDFDHCLSPNNHQCNLPKKLFNCDNFFTQCSKLTIQQCNYIQDNHYKSIDCDTMFNKTLNVTRRKIQHSKKSYKCNECGKAFKYCSSYRKHSIIHTGEKPYKCKLCGKSFTQCASLKKHQRIHTGEKPYRCEECGRSFNHYSILGQHQRIHTGEKPYKCKQCGKSFTQCSSLQKHQVIHTGEKPYRCAECGKSFTQNSTLSQHQRIHTGEKPYKCEECGKAFTQCSSLRKHQRIHTGEKPYKCEVCGRAFNCRSSFTKHKRIHTGEKPYKCKDCDKAFIHCTNLIQHQRIHTGEKPYKCNECGKSFSQCSNLRKHERIHT"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "17577",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "17346"
                }
            }
        ],
        "INSDFeature_key": "exon",
        "INSDFeature_location": "complement(17346..17577)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Zfp386"
                    },
                    {
                        "INSDQualifier_name": "number",
                        "INSDQualifier_value": "2"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "24522",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "24433"
                }
            }
        ],
        "INSDFeature_key": "exon",
        "INSDFeature_location": "complement(24433..24522)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Zfp386"
                    },
                    {
                        "INSDQualifier_name": "number",
                        "INSDQualifier_value": "1"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "62205",
                        "INSDInterval_to": "62250"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "62333",
                        "INSDInterval_to": "62639"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(62205..62250,62333..62639)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.89pg.195"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "62640",
                    "INSDInterval_to": "62646"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "62640..62646",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "62669",
                    "INSDInterval_to": "62677"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "62669..62677",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "71837",
                        "INSDInterval_to": "71882"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "71966",
                        "INSDInterval_to": "72270"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(71837..71882,71966..72270)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.88.194"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "72271",
                    "INSDInterval_to": "72277"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "72271..72277",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "72301",
                    "INSDInterval_to": "72309"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "72301..72309",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "78928",
                    "INSDInterval_to": "90035"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "78928..90035",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Odc1-like"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "91164",
                        "INSDInterval_to": "91209"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "91292",
                        "INSDInterval_to": "91596"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(91164..91209,91292..91596)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.87.193"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "91597",
                    "INSDInterval_to": "91603"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "91597..91603",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "91627",
                    "INSDInterval_to": "91635"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "91627..91635",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "108048",
                        "INSDInterval_to": "108093"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "108191",
                        "INSDInterval_to": "108481"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(108048..108093,108191..108481)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.86.192"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "108482",
                    "INSDInterval_to": "108488"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "108482..108488",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "108512",
                    "INSDInterval_to": "108520"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "108512..108520",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "119334",
                        "INSDInterval_to": "119379"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "119462",
                        "INSDInterval_to": "119766"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(119334..119379,119462..119766)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.85.191"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "119767",
                    "INSDInterval_to": "119773"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "119767..119773",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "119796",
                    "INSDInterval_to": "119804"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "119796..119804",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "151583",
                        "INSDInterval_to": "151628"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "151712",
                        "INSDInterval_to": "152016"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(151583..151628,151712..152016)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.84.190"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "152017",
                    "INSDInterval_to": "152023"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "152017..152023",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "152047",
                    "INSDInterval_to": "152055"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "152047..152055",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "159526",
                        "INSDInterval_to": "159566"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "159649",
                        "INSDInterval_to": "159953"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(159526..159566,159649..159953)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.83.189"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "159954",
                    "INSDInterval_to": "159960"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "159954..159960",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "159984",
                    "INSDInterval_to": "159992"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "159984..159992",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "179493",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "175243"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "complement(175243..179493)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Attp"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "183524",
                        "INSDInterval_to": "183569"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "183652",
                        "INSDInterval_to": "183956"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(183524..183569,183652..183956)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.82pg.188"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "183957",
                    "INSDInterval_to": "183963"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "183957..183963",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "183987",
                    "INSDInterval_to": "183995"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "183987..183995",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "202850",
                        "INSDInterval_to": "202895"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "202975",
                        "INSDInterval_to": "203279"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(202850..202895,202975..203279)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.81.187"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "203280",
                    "INSDInterval_to": "203286"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "203280..203286",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "203310",
                    "INSDInterval_to": "203318"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "203310..203318",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "209752",
                        "INSDInterval_to": "209797"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "209881",
                        "INSDInterval_to": "210185"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(209752..209797,209881..210185)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.80.186"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "210186",
                    "INSDInterval_to": "210192"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "210186..210192",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "210216",
                    "INSDInterval_to": "210224"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "210216..210224",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "213618",
                        "INSDInterval_to": "213663"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "213742",
                        "INSDInterval_to": "214051"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(213618..213663,213742..214051)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.16pg.185"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "223738",
                        "INSDInterval_to": "223783"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "223867",
                        "INSDInterval_to": "224171"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(223738..223783,223867..224171)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.79.184"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "224172",
                    "INSDInterval_to": "224178"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "224172..224178",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "224202",
                    "INSDInterval_to": "224210"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "224202..224210",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "227611",
                        "INSDInterval_to": "227656"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "227735",
                        "INSDInterval_to": "228044"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(227611..227656,227735..228044)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.15pg.183"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "237658",
                        "INSDInterval_to": "237703"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "237800",
                        "INSDInterval_to": "238104"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(237658..237703,237800..238104)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.78.182"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "238105",
                    "INSDInterval_to": "238111"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "238105..238111",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "238135",
                    "INSDInterval_to": "238143"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "238135..238143",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "240794",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "240067"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "complement(240067..240794)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Rpl7a-like"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "412141",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "256335"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "complement(256335..412141)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Hmgb2"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "263195",
                        "INSDInterval_to": "263240"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "263317",
                        "INSDInterval_to": "263628"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(263195..263240,263317..263628)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.14pg.181"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "268975",
                        "INSDInterval_to": "269020"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "269103",
                        "INSDInterval_to": "269407"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(268975..269020,269103..269407)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.77.180"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "269408",
                    "INSDInterval_to": "269414"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "269408..269414",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "269438",
                    "INSDInterval_to": "269446"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "269438..269446",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "281030",
                        "INSDInterval_to": "281075"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "281159",
                        "INSDInterval_to": "281464"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(281030..281075,281159..281464)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.76pg.179"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "281465",
                    "INSDInterval_to": "281471"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "281465..281471",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "281496",
                    "INSDInterval_to": "281503"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "281496..281503",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "305432",
                    "INSDInterval_to": "305743"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "305432..305743",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.13pg.178"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "312662",
                        "INSDInterval_to": "312707"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "312790",
                        "INSDInterval_to": "313094"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(312662..312707,312790..313094)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.75.177"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "313095",
                    "INSDInterval_to": "313101"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "313095..313101",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "313125",
                    "INSDInterval_to": "313133"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "313125..313133",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "328431",
                        "INSDInterval_to": "328476"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "328561",
                        "INSDInterval_to": "328873"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(328431..328476,328561..328873)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.74.176"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "328874",
                    "INSDInterval_to": "328880"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "328874..328880",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "328904",
                    "INSDInterval_to": "328912"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "328904..328912",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "378756",
                        "INSDInterval_to": "378801"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "378885",
                        "INSDInterval_to": "379189"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(378756..378801,378885..379189)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.73pg.175"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "379190",
                    "INSDInterval_to": "379196"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "379190..379196",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "379220",
                    "INSDInterval_to": "379228"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "379220..379228",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "422691",
                        "INSDInterval_to": "422736"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "422822",
                        "INSDInterval_to": "423133"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(422691..422736,422822..423133)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.12.174"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "423134",
                    "INSDInterval_to": "423140"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "423134..423140",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "423164",
                    "INSDInterval_to": "423172"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "423164..423172",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "447483",
                        "INSDInterval_to": "447528"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "447613",
                        "INSDInterval_to": "447917"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(447483..447528,447613..447917)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.72.173"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "447918",
                    "INSDInterval_to": "447924"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "447918..447924",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "447948",
                    "INSDInterval_to": "447956"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "447948..447956",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "459315",
                        "INSDInterval_to": "459360"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "459444",
                        "INSDInterval_to": "459737"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(459315..459360,459444..459737)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.71pg.172"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "466834",
                    "INSDInterval_to": "467138"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "466834..467138",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.70pg.171"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "467139",
                    "INSDInterval_to": "467145"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "467139..467145",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "467169",
                    "INSDInterval_to": "467177"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "467169..467177",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "477535",
                        "INSDInterval_to": "477580"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "477664",
                        "INSDInterval_to": "477968"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(477535..477580,477664..477968)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.69.170"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "477969",
                    "INSDInterval_to": "477975"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "477969..477975",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "477999",
                    "INSDInterval_to": "478007"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "477999..478007",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "480709",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "479805"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "complement(479805..480709)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Rpl7a-like"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "503487",
                        "INSDInterval_to": "503532"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "503618",
                        "INSDInterval_to": "503929"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(503487..503532,503618..503929)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.11.169"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "503930",
                    "INSDInterval_to": "503936"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "503930..503936",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "503960",
                    "INSDInterval_to": "503968"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "503960..503968",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "524683",
                    "INSDInterval_to": "525060"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "524683..525060",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Tceb2"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "538593",
                        "INSDInterval_to": "538638"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "538720",
                        "INSDInterval_to": "539022"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(538593..538638,538720..539022)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.68pg.168"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "539023",
                    "INSDInterval_to": "539029"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "539023..539029",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer - non consensus"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "539053",
                    "INSDInterval_to": "539061"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "539053..539061",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "897234",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "543980"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "complement(543980..897234)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Hmgb2"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "551296",
                        "INSDInterval_to": "551341"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "551427",
                        "INSDInterval_to": "551738"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(551296..551341,551427..551738)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.10pg.167"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "563101",
                        "INSDInterval_to": "563146"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "563229",
                        "INSDInterval_to": "563533"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(563101..563146,563229..563533)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.67.166"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "563534",
                    "INSDInterval_to": "563540"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "563534..563540",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "563564",
                    "INSDInterval_to": "563572"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "563564..563572",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "575012",
                        "INSDInterval_to": "575057"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "575141",
                        "INSDInterval_to": "575445"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(575012..575057,575141..575445)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.66.165"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "575446",
                    "INSDInterval_to": "575452"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "575446..575452",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "575476",
                    "INSDInterval_to": "575484"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "575476..575484",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "602296",
                        "INSDInterval_to": "602341"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "602428",
                        "INSDInterval_to": "602739"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(602296..602341,602428..602739)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.9.164"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "602740",
                    "INSDInterval_to": "602746"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "602740..602746",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "602770",
                    "INSDInterval_to": "602778"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "602770..602778",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "609640",
                        "INSDInterval_to": "609685"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "609768",
                        "INSDInterval_to": "610072"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(609640..609685,609768..610072)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.65.163"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "610073",
                    "INSDInterval_to": "610079"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "610073..610079",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "610103",
                    "INSDInterval_to": "610111"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "610103..610111",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "624221",
                        "INSDInterval_to": "624266"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "624350",
                        "INSDInterval_to": "624662"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(624221..624266,624350..624662)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.64.162"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "624663",
                    "INSDInterval_to": "624669"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "624663..624669",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "624693",
                    "INSDInterval_to": "624701"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "624693..624701",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "657956",
                    "INSDInterval_to": "658151"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "657956..658151",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.19.161"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "666359",
                    "INSDInterval_to": "666661"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "666359..666661",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.8pg.160"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "666662",
                    "INSDInterval_to": "666668"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "666662..666668",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "666692",
                    "INSDInterval_to": "666700"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "666692..666700",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "683950",
                        "INSDInterval_to": "683995"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "684079",
                        "INSDInterval_to": "684383"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(683950..683995,684079..684383)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.63pg.159"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "684384",
                    "INSDInterval_to": "684389"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "684384..684389",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "684413",
                    "INSDInterval_to": "684421"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "684413..684421",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "698645",
                        "INSDInterval_to": "698690"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "698789",
                        "INSDInterval_to": "699092"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(698645..698690,698789..699092)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.62pg.158"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "699093",
                    "INSDInterval_to": "699099"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "699093..699099",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "699123",
                    "INSDInterval_to": "699131"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "699123..699131",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "711506",
                        "INSDInterval_to": "711551"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "711634",
                        "INSDInterval_to": "711938"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(711506..711551,711634..711938)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.61.157"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "711939",
                    "INSDInterval_to": "711945"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "711939..711945",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "711969",
                    "INSDInterval_to": "711977"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "711969..711977",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "726115",
                        "INSDInterval_to": "726157"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "726237",
                        "INSDInterval_to": "726547"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(726115..726157,726237..726547)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.60pg.156"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "726548",
                    "INSDInterval_to": "726554"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "726548..726554",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "726578",
                    "INSDInterval_to": "726586"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "726578..726586",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "735564",
                        "INSDInterval_to": "735609"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "735692",
                        "INSDInterval_to": "735996"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(735564..735609,735692..735996)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.59.155"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "735997",
                    "INSDInterval_to": "736003"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "735997..736003",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "736027",
                    "INSDInterval_to": "736035"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "736027..736035",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "758479",
                        "INSDInterval_to": "758524"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "758608",
                        "INSDInterval_to": "758912"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(758479..758524,758608..758912)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.58.154"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "758913",
                    "INSDInterval_to": "758919"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "758913..758919",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "758943",
                    "INSDInterval_to": "758951"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "758943..758951",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "776574",
                        "INSDInterval_to": "776619"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "776705",
                        "INSDInterval_to": "777016"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(776574..776619,776705..777016)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.7.153"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "777017",
                    "INSDInterval_to": "777023"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "777017..777023",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "777047",
                    "INSDInterval_to": "777055"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "777047..777055",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "802222",
                        "INSDInterval_to": "802267"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "802352",
                        "INSDInterval_to": "802654"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(802222..802267,802352..802654)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.57pg.152"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "802655",
                    "INSDInterval_to": "802661"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "802655..802661",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "802685",
                    "INSDInterval_to": "802693"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "802685..802693",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "812649",
                    "INSDInterval_to": "812941"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "812649..812941",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.6pg.151"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "812947",
                    "INSDInterval_to": "812953"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "812947..812953",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "812977",
                    "INSDInterval_to": "812985"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "812977..812985",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "826533",
                    "INSDInterval_to": "826895"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "826533..826895",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Birc3"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "827843",
                        "INSDInterval_to": "827888"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "827972",
                        "INSDInterval_to": "828276"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(827843..827888,827972..828276)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.56.150"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "828277",
                    "INSDInterval_to": "828283"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "828277..828283",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "828307",
                    "INSDInterval_to": "828315"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "828307..828315",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "862551",
                        "INSDInterval_to": "862596"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "862679",
                        "INSDInterval_to": "862983"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(862551..862596,862679..862983)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.55.149"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "862984",
                    "INSDInterval_to": "862990"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "862984..862990",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "863014",
                    "INSDInterval_to": "863022"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "863014..863022",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "876975",
                        "INSDInterval_to": "877015"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "877110",
                        "INSDInterval_to": "877403"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(876975..877015,877110..877403)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.54.148"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "877404",
                    "INSDInterval_to": "877410"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "877404..877410",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "877434",
                    "INSDInterval_to": "877442"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "877434..877442",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "904858",
                        "INSDInterval_to": "904903"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "904990",
                        "INSDInterval_to": "905301"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(904858..904903,904990..905301)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.5.147"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "905302",
                    "INSDInterval_to": "905308"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "905302..905308",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "905332",
                    "INSDInterval_to": "905340"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "905332..905340",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "912243",
                        "INSDInterval_to": "912288"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "912371",
                        "INSDInterval_to": "912675"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(912243..912288,912371..912675)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.53.146"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "912676",
                    "INSDInterval_to": "912682"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "912676..912682",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "912706",
                    "INSDInterval_to": "912714"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "912706..912714",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "925159",
                        "INSDInterval_to": "925204"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "925287",
                        "INSDInterval_to": "925591"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(925159..925204,925287..925591)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.52.145"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "925592",
                    "INSDInterval_to": "925598"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "925592..925598",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "925622",
                    "INSDInterval_to": "925630"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "925622..925630",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "939277",
                        "INSDInterval_to": "939319"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "939399",
                        "INSDInterval_to": "939709"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(939277..939319,939399..939709)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.51pg.144"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "939710",
                    "INSDInterval_to": "939716"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "939710..939716",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "939740",
                    "INSDInterval_to": "939748"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "939740..939748",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "950898",
                        "INSDInterval_to": "950943"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "951026",
                        "INSDInterval_to": "951330"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(950898..950943,951026..951330)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.50.143"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "951331",
                    "INSDInterval_to": "951337"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "951331..951337",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "951361",
                    "INSDInterval_to": "951369"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "951361..951369",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1003076",
                        "INSDInterval_to": "1003121"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1003207",
                        "INSDInterval_to": "1003518"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1003076..1003121,1003207..1003518)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.4.142"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1003519",
                    "INSDInterval_to": "1003525"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1003519..1003525",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1003549",
                    "INSDInterval_to": "1003557"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1003549..1003557",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1015421",
                        "INSDInterval_to": "1015466"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1015550",
                        "INSDInterval_to": "1015854"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1015421..1015466,1015550..1015854)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.49.141"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1015855",
                    "INSDInterval_to": "1015861"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1015855..1015861",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1015885",
                    "INSDInterval_to": "1015893"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1015885..1015893",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1032423",
                        "INSDInterval_to": "1032461"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1032546",
                        "INSDInterval_to": "1032850"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1032423..1032461,1032546..1032850)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.48pg.140"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1032851",
                    "INSDInterval_to": "1032857"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1032851..1032857",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1032881",
                    "INSDInterval_to": "1032889"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1032881..1032889",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1046618",
                        "INSDInterval_to": "1046663"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1046747",
                        "INSDInterval_to": "1047051"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1046618..1046663,1046747..1047051)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.3.139"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1047052",
                    "INSDInterval_to": "1047058"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1047052..1047058",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1047082",
                    "INSDInterval_to": "1047090"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1047082..1047090",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1067380",
                    "INSDInterval_to": "1067651"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "1067380..1067651",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.2pg.138"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1079543",
                        "INSDInterval_to": "1079588"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1079665",
                        "INSDInterval_to": "1079969"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1079543..1079588,1079665..1079969)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.47.137"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1079970",
                    "INSDInterval_to": "1079976"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1079970..1079976",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1080000",
                    "INSDInterval_to": "1080008"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1080000..1080008",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1092951",
                        "INSDInterval_to": "1092996"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1093080",
                        "INSDInterval_to": "1093384"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1092951..1092996,1093080..1093384)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.46pg.136"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1093385",
                    "INSDInterval_to": "1093391"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1093385..1093391",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1093415",
                    "INSDInterval_to": "1093421"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1093415..1093421",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1100717",
                    "INSDInterval_to": "1101022"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "1100717..1101022",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.45pg.135"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1101023",
                    "INSDInterval_to": "1101029"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1101023..1101029",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1101053",
                    "INSDInterval_to": "1101061"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1101053..1101061",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1107264",
                        "INSDInterval_to": "1107309"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1107392",
                        "INSDInterval_to": "1107697"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1107264..1107309,1107392..1107697)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.44pg.134"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1107698",
                    "INSDInterval_to": "1107704"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1107698..1107704",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1107728",
                    "INSDInterval_to": "1107736"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1107728..1107736",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1124694",
                        "INSDInterval_to": "1124739"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1124823",
                        "INSDInterval_to": "1125127"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1124694..1124739,1124823..1125127)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.43.133"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1125128",
                    "INSDInterval_to": "1125134"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1125128..1125134",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1125158",
                    "INSDInterval_to": "1125166"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1125158..1125166",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1133533",
                        "INSDInterval_to": "1133578"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1133660",
                        "INSDInterval_to": "1133964"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1133533..1133578,1133660..1133964)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.42.132"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1133965",
                    "INSDInterval_to": "1133971"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1133965..1133971",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1133995",
                    "INSDInterval_to": "1134003"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1133995..1134003",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1138295",
                        "INSDInterval_to": "1138340"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1138424",
                        "INSDInterval_to": "1138706"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1138295..1138340,1138424..1138706)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.41pg.131"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1138707",
                    "INSDInterval_to": "1138713"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1138707..1138713",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1138737",
                    "INSDInterval_to": "1138745"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1138737..1138745",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1145704",
                        "INSDInterval_to": "1145749"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1145833",
                        "INSDInterval_to": "1146136"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1145704..1145749,1145833..1146136)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.40pg.130"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1146137",
                    "INSDInterval_to": "1146143"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1146137..1146143",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1146167",
                    "INSDInterval_to": "1146175"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1146167..1146175",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1156045",
                        "INSDInterval_to": "1156090"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1156174",
                        "INSDInterval_to": "1156478"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1156045..1156090,1156174..1156478)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.39.129"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1156479",
                    "INSDInterval_to": "1156485"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1156479..1156485",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1156509",
                    "INSDInterval_to": "1156517"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1156509..1156517",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1160634",
                        "INSDInterval_to": "1160679"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1160763",
                        "INSDInterval_to": "1161066"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1160634..1160679,1160763..1161066)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.38pg.128"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1161067",
                    "INSDInterval_to": "1161073"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1161067..1161073",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1161096",
                    "INSDInterval_to": "1161104"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1161096..1161104",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1174406",
                        "INSDInterval_to": "1174451"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1174535",
                        "INSDInterval_to": "1174839"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1174406..1174451,1174535..1174839)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.37.127"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1174840",
                    "INSDInterval_to": "1174846"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1174840..1174846",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1174869",
                    "INSDInterval_to": "1174877"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1174869..1174877",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1190756",
                        "INSDInterval_to": "1190801"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1190885",
                        "INSDInterval_to": "1191189"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1190756..1190801,1190885..1191189)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.36.126"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1191190",
                    "INSDInterval_to": "1191196"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1191190..1191196",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1191219",
                    "INSDInterval_to": "1191227"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1191219..1191227",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1195341",
                        "INSDInterval_to": "1195386"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1195470",
                        "INSDInterval_to": "1195774"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1195341..1195386,1195470..1195774)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.35pg.125"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1195775",
                    "INSDInterval_to": "1195781"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1195775..1195781",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer - non consensus"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1195805",
                    "INSDInterval_to": "1195813"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1195805..1195813",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1219457",
                        "INSDInterval_to": "1219502"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1219586",
                        "INSDInterval_to": "1219890"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1219457..1219502,1219586..1219890)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.34.124"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1219891",
                    "INSDInterval_to": "1219897"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1219891..1219897",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1219920",
                    "INSDInterval_to": "1219928"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1219920..1219928",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1226997",
                        "INSDInterval_to": "1227042"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1227126",
                        "INSDInterval_to": "1227429"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1226997..1227042,1227126..1227429)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.33pg.123"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1227430",
                    "INSDInterval_to": "1227436"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1227430..1227436",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1227459",
                    "INSDInterval_to": "1227467"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1227459..1227467",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1235424",
                        "INSDInterval_to": "1235469"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1235553",
                        "INSDInterval_to": "1235857"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1235424..1235469,1235553..1235857)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.32pg.122"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1235858",
                    "INSDInterval_to": "1235864"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1235858..1235864",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1235887",
                    "INSDInterval_to": "1235896"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1235887..1235896",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1241383",
                        "INSDInterval_to": "1241428"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1241512",
                        "INSDInterval_to": "1241816"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1241383..1241428,1241512..1241816)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.31.121"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1241817",
                    "INSDInterval_to": "1241823"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1241817..1241823",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1241847",
                    "INSDInterval_to": "1241855"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1241847..1241855",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1253554",
                        "INSDInterval_to": "1253599"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1253679",
                        "INSDInterval_to": "1253980"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1253554..1253599,1253679..1253980)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.30pg.120"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1253982",
                    "INSDInterval_to": "1253988"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1253982..1253988",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1254009",
                    "INSDInterval_to": "1254017"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1254009..1254017",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1258331",
                        "INSDInterval_to": "1258380"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1258464",
                        "INSDInterval_to": "1258759"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1258331..1258380,1258464..1258759)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.29pg.119"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1258760",
                    "INSDInterval_to": "1258766"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1258760..1258766",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1258790",
                    "INSDInterval_to": "1258798"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1258790..1258798",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1261226",
                        "INSDInterval_to": "1261271"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1261355",
                        "INSDInterval_to": "1261659"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1261226..1261271,1261355..1261659)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.28pg.118"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1261660",
                    "INSDInterval_to": "1261666"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1261660..1261666",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1261689",
                    "INSDInterval_to": "1261697"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1261689..1261697",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1266032",
                        "INSDInterval_to": "1266077"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1266160",
                        "INSDInterval_to": "1266464"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1266032..1266077,1266160..1266464)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.27pg.117"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1266465",
                    "INSDInterval_to": "1266471"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1266465..1266471",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1266495",
                    "INSDInterval_to": "1266503"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1266495..1266503",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1282275",
                        "INSDInterval_to": "1282320"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1282404",
                        "INSDInterval_to": "1282708"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1282275..1282320,1282404..1282708)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.26.116"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1282709",
                    "INSDInterval_to": "1282715"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1282709..1282715",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1282739",
                    "INSDInterval_to": "1282747"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1282739..1282747",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1289245",
                        "INSDInterval_to": "1289285"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1289367",
                        "INSDInterval_to": "1289671"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1289245..1289285,1289367..1289671)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.25pg.115"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1289672",
                    "INSDInterval_to": "1289678"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1289672..1289678",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1289702",
                    "INSDInterval_to": "1289710"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1289702..1289710",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1297719",
                        "INSDInterval_to": "1297764"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1297848",
                        "INSDInterval_to": "1298152"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1297719..1297764,1297848..1298152)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.24pg.114"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1298180",
                    "INSDInterval_to": "1298188"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1298180..1298188",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer - no heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1306197",
                        "INSDInterval_to": "1306242"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1306326",
                        "INSDInterval_to": "1306630"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1306197..1306242,1306326..1306630)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.23.113"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1306631",
                    "INSDInterval_to": "1306637"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1306631..1306637",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1306661",
                    "INSDInterval_to": "1306669"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1306661..1306669",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1324374",
                        "INSDInterval_to": "1324419"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1324503",
                        "INSDInterval_to": "1324807"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1324374..1324419,1324503..1324807)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.22.112"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1324808",
                    "INSDInterval_to": "1324814"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1324808..1324814",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1324838",
                    "INSDInterval_to": "1324846"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1324838..1324846",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1330244",
                        "INSDInterval_to": "1330284"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1330368",
                        "INSDInterval_to": "1330672"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1330244..1330284,1330368..1330672)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.21pg.111"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1330673",
                    "INSDInterval_to": "1330679"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1330673..1330679",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1330703",
                    "INSDInterval_to": "1330711"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1330703..1330711",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1334402",
                        "INSDInterval_to": "1334447"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1334531",
                        "INSDInterval_to": "1334835"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1334402..1334447,1334531..1334835)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.20pg.110"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1334836",
                    "INSDInterval_to": "1334842"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1334836..1334842",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1334865",
                    "INSDInterval_to": "1334873"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1334865..1334873",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1346875",
                        "INSDInterval_to": "1346920"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1347004",
                        "INSDInterval_to": "1347308"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1346875..1346920,1347004..1347308)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.19.109"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1347309",
                    "INSDInterval_to": "1347315"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1347309..1347315",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1347339",
                    "INSDInterval_to": "1347347"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1347339..1347347",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1361999",
                        "INSDInterval_to": "1362044"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1362128",
                        "INSDInterval_to": "1362432"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1361999..1362044,1362128..1362432)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.18.108"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1362433",
                    "INSDInterval_to": "1362439"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1362433..1362439",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1362462",
                    "INSDInterval_to": "1362470"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1362462..1362470",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1366500",
                        "INSDInterval_to": "1366545"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1366629",
                        "INSDInterval_to": "1366933"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1366500..1366545,1366629..1366933)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.17pg.107"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1366934",
                    "INSDInterval_to": "1366940"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1366934..1366940",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1366964",
                    "INSDInterval_to": "1366972"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1366964..1366972",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1388043",
                        "INSDInterval_to": "1388088"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1388172",
                        "INSDInterval_to": "1388476"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1388043..1388088,1388172..1388476)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.16.106"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1388477",
                    "INSDInterval_to": "1388483"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1388477..1388483",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1388507",
                    "INSDInterval_to": "1388515"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1388507..1388515",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1393927",
                        "INSDInterval_to": "1393972"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1394056",
                        "INSDInterval_to": "1394359"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1393927..1393972,1394056..1394359)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.15pg.105"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1394360",
                    "INSDInterval_to": "1394366"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1394360..1394366",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1394390",
                    "INSDInterval_to": "1394398"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1394390..1394398",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1398707",
                        "INSDInterval_to": "1398753"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1398836",
                        "INSDInterval_to": "1399140"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1398707..1398753,1398836..1399140)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.14pg.104"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1399141",
                    "INSDInterval_to": "1399147"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1399141..1399147",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1399170",
                    "INSDInterval_to": "1399178"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1399170..1399178",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1404808",
                        "INSDInterval_to": "1404848"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1404932",
                        "INSDInterval_to": "1405236"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1404808..1404848,1404932..1405236)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.13.103"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1405237",
                    "INSDInterval_to": "1405243"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1405237..1405243",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1405266",
                    "INSDInterval_to": "1405274"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1405266..1405274",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1413322",
                        "INSDInterval_to": "1413367"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1413451",
                        "INSDInterval_to": "1413755"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1413322..1413367,1413451..1413755)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.12.102"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1413756",
                    "INSDInterval_to": "1413762"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1413756..1413762",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1413786",
                    "INSDInterval_to": "1413794"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1413786..1413794",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1424103",
                        "INSDInterval_to": "1424148"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1424232",
                        "INSDInterval_to": "1424536"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1424103..1424148,1424232..1424536)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.11pg.101"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1424537",
                    "INSDInterval_to": "1424543"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1424537..1424543",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1424567",
                    "INSDInterval_to": "1424575"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1424567..1424575",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1439997",
                        "INSDInterval_to": "1440041"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1440138",
                        "INSDInterval_to": "1440428"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1439997..1440041,1440138..1440428)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.10pg.100"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1440429",
                    "INSDInterval_to": "1440435"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1440429..1440435",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1440459",
                    "INSDInterval_to": "1440467"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1440459..1440467",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1454825",
                        "INSDInterval_to": "1454870"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1454954",
                        "INSDInterval_to": "1455258"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1454825..1454870,1454954..1455258)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.9.99"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1455259",
                    "INSDInterval_to": "1455265"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1455259..1455265",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1455289",
                    "INSDInterval_to": "1455297"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1455289..1455297",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1458433",
                        "INSDInterval_to": "1458478"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1458561",
                        "INSDInterval_to": "1458865"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1458433..1458478,1458561..1458865)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.8.98"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1458866",
                    "INSDInterval_to": "1458872"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1458866..1458872",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1458896",
                    "INSDInterval_to": "1458904"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1458896..1458904",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1472782",
                        "INSDInterval_to": "1472827"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1472911",
                        "INSDInterval_to": "1473228"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1472782..1472827,1472911..1473228)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.7pg.97"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1473230",
                    "INSDInterval_to": "1473236"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1473230..1473236",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1473260",
                    "INSDInterval_to": "1473268"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1473260..1473268",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1487106",
                        "INSDInterval_to": "1487151"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1487235",
                        "INSDInterval_to": "1487539"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1487106..1487151,1487235..1487539)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.6.96"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1487540",
                    "INSDInterval_to": "1487546"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1487540..1487546",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1487570",
                    "INSDInterval_to": "1487578"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1487570..1487578",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1506076",
                        "INSDInterval_to": "1506145"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1506223",
                        "INSDInterval_to": "1506530"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1506076..1506145,1506223..1506530)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VH15.1.95"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1506531",
                    "INSDInterval_to": "1506537"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1506531..1506537",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1506561",
                    "INSDInterval_to": "1506569"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1506561..1506569",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1515065",
                        "INSDInterval_to": "1515110"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1515195",
                        "INSDInterval_to": "1515499"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1515065..1515110,1515195..1515499)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.5pg.94"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1515500",
                    "INSDInterval_to": "1515506"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1515500..1515506",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1515530",
                    "INSDInterval_to": "1515538"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1515530..1515538",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1532180",
                        "INSDInterval_to": "1532225"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1532309",
                        "INSDInterval_to": "1532613"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1532180..1532225,1532309..1532613)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.4.93"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1532614",
                    "INSDInterval_to": "1532620"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1532614..1532620",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1532644",
                    "INSDInterval_to": "1532652"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1532644..1532652",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1537848",
                    "INSDInterval_to": "1538091"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "1537848..1538091",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.18.92"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1547222",
                        "INSDInterval_to": "1547267"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1547356",
                        "INSDInterval_to": "1547668"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1547222..1547267,1547356..1547668)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VH10.3.91"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1547669",
                    "INSDInterval_to": "1547675"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1547669..1547675",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1547699",
                    "INSDInterval_to": "1547707"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1547699..1547707",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1557346",
                        "INSDInterval_to": "1557391"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1557475",
                        "INSDInterval_to": "1557779"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1557346..1557391,1557475..1557779)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.3.90"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1557780",
                    "INSDInterval_to": "1557786"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1557780..1557786",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1557810",
                    "INSDInterval_to": "1557818"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1557810..1557818",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1574806",
                    "INSDInterval_to": "1575106"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "1574806..1575106",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VH10.2pg.89"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1575107",
                    "INSDInterval_to": "1575113"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1575107..1575113",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer - non consensus"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1575137",
                    "INSDInterval_to": "1575145"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1575137..1575145",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1583540",
                        "INSDInterval_to": "1583585"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1583669",
                        "INSDInterval_to": "1583973"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1583540..1583585,1583669..1583973)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.2.88"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1583974",
                    "INSDInterval_to": "1583980"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1583974..1583980",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1584003",
                    "INSDInterval_to": "1584011"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1584003..1584011",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1589192",
                    "INSDInterval_to": "1589438"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "1589192..1589438",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.17.87"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1591658",
                        "INSDInterval_to": "1591703"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1591792",
                        "INSDInterval_to": "1592104"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1591658..1591703,1591792..1592104)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VH10.1.86"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1592105",
                    "INSDInterval_to": "1592111"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1592105..1592111",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1592135",
                    "INSDInterval_to": "1592143"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1592135..1592143",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1601574",
                        "INSDInterval_to": "1601619"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1601703",
                        "INSDInterval_to": "1602034"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1601574..1601619,1601703..1602034)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J558.1.85"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1602035",
                    "INSDInterval_to": "1602041"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1602035..1602041",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1602065",
                    "INSDInterval_to": "1602073"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1602065..1602073",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1608379",
                        "INSDInterval_to": "1608427"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1608523",
                        "INSDInterval_to": "1608819"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1608379..1608427,1608523..1608819)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609.1.84"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1608820",
                    "INSDInterval_to": "1608826"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1608820..1608826",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1608849",
                    "INSDInterval_to": "1608857"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1608849..1608857",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1615026",
                        "INSDInterval_to": "1615071"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1615173",
                        "INSDInterval_to": "1615483"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1615026..1615071,1615173..1615483)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J606.5.83"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1615484",
                    "INSDInterval_to": "1615490"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1615484..1615490",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1615514",
                    "INSDInterval_to": "1615522"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1615514..1615522",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1635863",
                        "INSDInterval_to": "1635908"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1636009",
                        "INSDInterval_to": "1636319"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1635863..1635908,1636009..1636319)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J606.4.82"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1636320",
                    "INSDInterval_to": "1636326"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1636320..1636326",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1636350",
                    "INSDInterval_to": "1636358"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1636350..1636358",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1654056",
                        "INSDInterval_to": "1654110"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1654203",
                        "INSDInterval_to": "1654513"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1654056..1654110,1654203..1654513)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J606.3.81"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1654514",
                    "INSDInterval_to": "1654520"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1654514..1654520",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1654544",
                    "INSDInterval_to": "1654552"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1654544..1654552",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1664148",
                        "INSDInterval_to": "1664193"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1664325",
                        "INSDInterval_to": "1664635"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1664148..1664193,1664325..1664635)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J606.2.80"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1664636",
                    "INSDInterval_to": "1664642"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1664636..1664642",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1664666",
                    "INSDInterval_to": "1664674"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1664666..1664674",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1678940",
                        "INSDInterval_to": "1678994"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1679088",
                        "INSDInterval_to": "1679397"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1678940..1678994,1679088..1679397)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "J606.1.79"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1679398",
                    "INSDInterval_to": "1679404"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1679398..1679404",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1679428",
                    "INSDInterval_to": "1679436"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1679428..1679436",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1704155",
                        "INSDInterval_to": "1704197"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1704279",
                        "INSDInterval_to": "1704589"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1704155..1704197,1704279..1704589)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VH12.1.78"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1704590",
                    "INSDInterval_to": "1704596"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1704590..1704596",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1704620",
                    "INSDInterval_to": "1704628"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1704620..1704628",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1712891",
                        "INSDInterval_to": "1712936"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1713038",
                        "INSDInterval_to": "1713348"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1712891..1712936,1713038..1713348)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609N.2.77"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1713349",
                    "INSDInterval_to": "1713355"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1713349..1713355",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1713379",
                    "INSDInterval_to": "1713387"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1713379..1713387",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1719103",
                    "INSDInterval_to": "1719406"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "1719103..1719406",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.16.76"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1738026",
                    "INSDInterval_to": "1738306"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "1738026..1738306",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.15.75"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1748308",
                        "INSDInterval_to": "1748350"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1748432",
                        "INSDInterval_to": "1748735"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1748308..1748350,1748432..1748735)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "36-60.8.74"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1748736",
                    "INSDInterval_to": "1748742"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1748736..1748742",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1748766",
                    "INSDInterval_to": "1748774"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1748766..1748774",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1750644",
                        "INSDInterval_to": "1750688"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1750782",
                        "INSDInterval_to": "1751078"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1750644..1750688,1750782..1751078)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.14.73"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1751081",
                    "INSDInterval_to": "1751087"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1751081..1751087",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1751111",
                    "INSDInterval_to": "1751119"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1751111..1751119",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1756491",
                        "INSDInterval_to": "1756533"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1756619",
                        "INSDInterval_to": "1756925"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1756491..1756533,1756619..1756925)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "36-60.7pg.72"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1756926",
                    "INSDInterval_to": "1756932"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1756926..1756932",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer - non consensus"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1756956",
                    "INSDInterval_to": "1756964"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1756956..1756964",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1770720",
                        "INSDInterval_to": "1770765"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1770844",
                        "INSDInterval_to": "1771148"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1770720..1770765,1770844..1771148)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VGAM3.8-4-71"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1771149",
                    "INSDInterval_to": "1771155"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1771149..1771155",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1771178",
                    "INSDInterval_to": "1771186"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1771178..1771186",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1782528",
                        "INSDInterval_to": "1782570"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1782652",
                        "INSDInterval_to": "1782958"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1782528..1782570,1782652..1782958)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "36-60.6.70"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1782959",
                    "INSDInterval_to": "1782965"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1782959..1782965",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1782989",
                    "INSDInterval_to": "1782997"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1782989..1782997",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1786119",
                    "INSDInterval_to": "1786429"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "1786119..1786429",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.13.69"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1803046",
                        "INSDInterval_to": "1803090"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1803195",
                        "INSDInterval_to": "1803505"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1803046..1803090,1803195..1803505)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "3609N.1pg.68"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1803510",
                    "INSDInterval_to": "1803516"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1803510..1803516",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1803544",
                    "INSDInterval_to": "1803552"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1803544..1803552",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer - RSS spacer 27bp"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1804775",
                    "INSDInterval_to": "1805505"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "1804775..1805505",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Rad23b"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1808025",
                        "INSDInterval_to": "1808067"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1808151",
                        "INSDInterval_to": "1808460"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1808025..1808067,1808151..1808460)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "36-60.5.67"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1808461",
                    "INSDInterval_to": "1808467"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1808461..1808467",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1808491",
                    "INSDInterval_to": "1808499"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1808491..1808499",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1817061",
                        "INSDInterval_to": "1817103"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1817186",
                        "INSDInterval_to": "1817495"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1817061..1817103,1817186..1817495)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "36-60.4.66"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1817496",
                    "INSDInterval_to": "1817502"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1817496..1817502",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1817526",
                    "INSDInterval_to": "1817534"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1817526..1817534",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1847858",
                        "INSDInterval_to": "1847903"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1848008",
                        "INSDInterval_to": "1848318"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1847858..1847903,1848008..1848318)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "S107.4.65"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1848325",
                    "INSDInterval_to": "1848331"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1848325..1848331",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1848355",
                    "INSDInterval_to": "1848363"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1848355..1848363",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1874237",
                        "INSDInterval_to": "1874279"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1874367",
                        "INSDInterval_to": "1874673"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1874237..1874279,1874367..1874673)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "36-60.3.64"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1874674",
                    "INSDInterval_to": "1874681"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1874674..1874681",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer - 1bp insertion"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1874704",
                    "INSDInterval_to": "1874712"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1874704..1874712",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1882119",
                    "INSDInterval_to": "1882424"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "1882119..1882424",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Wdr22"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1894245",
                        "INSDInterval_to": "1894290"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1894384",
                        "INSDInterval_to": "1894674"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1894245..1894290,1894384..1894674)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "SM7.4.63"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1894675",
                    "INSDInterval_to": "1894681"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1894675..1894681",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1894704",
                    "INSDInterval_to": "1894712"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1894704..1894712",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1917468",
                        "INSDInterval_to": "1917513"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1917618",
                        "INSDInterval_to": "1917928"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1917468..1917513,1917618..1917928)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "S107.3.62"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1917933",
                    "INSDInterval_to": "1917939"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1917933..1917939",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1917963",
                    "INSDInterval_to": "1917971"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1917963..1917971",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1929990",
                        "INSDInterval_to": "1930035"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1930116",
                        "INSDInterval_to": "1930420"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1929990..1930035,1930116..1930420)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VGAM3.8-3-61"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1930421",
                    "INSDInterval_to": "1930427"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1930421..1930427",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1930450",
                    "INSDInterval_to": "1930458"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1930450..1930458",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1943266",
                    "INSDInterval_to": "1943571"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "1943266..1943571",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.12.60"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1961682",
                        "INSDInterval_to": "1961727"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1961807",
                        "INSDInterval_to": "1962111"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1961682..1961727,1961807..1962111)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VGAM3.8-2-59"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1962112",
                    "INSDInterval_to": "1962118"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1962112..1962118",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1962141",
                    "INSDInterval_to": "1962149"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1962141..1962149",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1963410",
                        "INSDInterval_to": "1963459"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1963538",
                        "INSDInterval_to": "1963853"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1963410..1963459,1963538..1963853)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.11.58"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1963854",
                    "INSDInterval_to": "1963860"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1963854..1963860",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1963884",
                    "INSDInterval_to": "1963892"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1963884..1963892",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1976754",
                        "INSDInterval_to": "1976799"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "1976880",
                        "INSDInterval_to": "1977184"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(1976754..1976799,1976880..1977184)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VGAM3.8-1-57"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1977185",
                    "INSDInterval_to": "1977191"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1977185..1977191",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1977214",
                    "INSDInterval_to": "1977222"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "1977214..1977222",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "1980658",
                    "INSDInterval_to": "1980945"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "1980658..1980945",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.10.56"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2001831",
                        "INSDInterval_to": "2001884"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2001975",
                        "INSDInterval_to": "2002284"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2001831..2001884,2001975..2002284)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VH16.1.55"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2002285",
                    "INSDInterval_to": "2002291"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2002285..2002291",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2002315",
                    "INSDInterval_to": "2002323"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2002315..2002323",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2010839",
                        "INSDInterval_to": "2010884"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2010963",
                        "INSDInterval_to": "2011267"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2010839..2010884,2010963..2011267)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "SM7.3.54"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2011268",
                    "INSDInterval_to": "2011274"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2011268..2011274",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2011297",
                    "INSDInterval_to": "2011305"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2011297..2011305",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2022408",
                        "INSDInterval_to": "2022453"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2022565",
                        "INSDInterval_to": "2022871"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2022408..2022453,2022565..2022871)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VH11.2.53"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2022872",
                    "INSDInterval_to": "2022878"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2022872..2022878",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2022902",
                    "INSDInterval_to": "2022910"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2022902..2022910",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2032867",
                    "INSDInterval_to": "2033174"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2032867..2033174",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.9.52"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2036879",
                        "INSDInterval_to": "2036921"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2037004",
                        "INSDInterval_to": "2037309"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2036879..2036921,2037004..2037309)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "36-60.2pg.51"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2037310",
                    "INSDInterval_to": "2037316"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2037310..2037316",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2037340",
                    "INSDInterval_to": "2037348"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2037340..2037348",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2057518",
                        "INSDInterval_to": "2057560"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2057662",
                        "INSDInterval_to": "2057966"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2057518..2057560,2057662..2057966)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "X24.2.50"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2057969",
                    "INSDInterval_to": "2057975"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2057969..2057975",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2057999",
                    "INSDInterval_to": "2058007"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2057999..2058007",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2076214",
                        "INSDInterval_to": "2076259"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2076339",
                        "INSDInterval_to": "2076643"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2076214..2076259,2076339..2076643)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "SM7.2.49"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2076644",
                    "INSDInterval_to": "2076650"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2076644..2076650",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2076673",
                    "INSDInterval_to": "2076681"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2076673..2076681",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2088771",
                        "INSDInterval_to": "2088816"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2088927",
                        "INSDInterval_to": "2089233"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2088771..2088816,2088927..2089233)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "VH11.1.48"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2089234",
                    "INSDInterval_to": "2089240"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2089234..2089240",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2089264",
                    "INSDInterval_to": "2089272"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2089264..2089272",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2102315",
                    "INSDInterval_to": "2102619"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2102315..2102619",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.8.47"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2106334",
                        "INSDInterval_to": "2106376"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2106458",
                        "INSDInterval_to": "2106764"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2106334..2106376,2106458..2106764)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "36-60.1.46"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2106765",
                    "INSDInterval_to": "2106771"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2106765..2106771",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2106795",
                    "INSDInterval_to": "2106803"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2106795..2106803",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2122417",
                        "INSDInterval_to": "2122459"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2122561",
                        "INSDInterval_to": "2122865"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2122417..2122459,2122561..2122865)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "X24.1pg.45"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2122868",
                    "INSDInterval_to": "2122874"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2122868..2122874",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2122898",
                    "INSDInterval_to": "2122906"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2122898..2122906",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2138767",
                        "INSDInterval_to": "2138812"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2138892",
                        "INSDInterval_to": "2139196"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2138767..2138812,2138892..2139196)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "SM7.1.44"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2139197",
                    "INSDInterval_to": "2139203"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2139197..2139203",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2139226",
                    "INSDInterval_to": "2139234"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2139226..2139234",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2158666",
                        "INSDInterval_to": "2158705"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2158808",
                        "INSDInterval_to": "2159116"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2158666..2158705,2158808..2159116)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "S107.2pg.43"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2159125",
                    "INSDInterval_to": "2159131"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2159125..2159131",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2159147",
                    "INSDInterval_to": "2159155"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2159147..2159155",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2174203",
                        "INSDInterval_to": "2174248"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2174425",
                        "INSDInterval_to": "2174739"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2174203..2174248,2174425..2174739)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "S107.1.42"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2174742",
                    "INSDInterval_to": "2174748"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2174742..2174748",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2174772",
                    "INSDInterval_to": "2174780"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2174772..2174780",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2186090",
                    "INSDInterval_to": "2186386"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2186090..2186386",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.7.41"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2191620",
                        "INSDInterval_to": "2191665"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2191750",
                        "INSDInterval_to": "2192053"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2191620..2191665,2191750..2192053)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.13.40"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2192054",
                    "INSDInterval_to": "2192060"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2192054..2192060",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2192084",
                    "INSDInterval_to": "2192092"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2192084..2192092",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2193498",
                    "INSDInterval_to": "2193794"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2193498..2193794",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.12pg.39"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2194910",
                        "INSDInterval_to": "2194955"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2195074",
                        "INSDInterval_to": "2195377"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2194910..2194955,2195074..2195377)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.21pg.38"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2211533",
                        "INSDInterval_to": "2211578"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2211696",
                        "INSDInterval_to": "2212000"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2211533..2211578,2211696..2212000)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.20.37"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2212001",
                    "INSDInterval_to": "2212007"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2212001..2212007",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2212031",
                    "INSDInterval_to": "2212039"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2212031..2212039",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2232166",
                        "INSDInterval_to": "2232211"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2232317",
                        "INSDInterval_to": "2232621"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2232166..2232211,2232317..2232621)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.19.36"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2232624",
                    "INSDInterval_to": "2232630"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2232624..2232630",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2232654",
                    "INSDInterval_to": "2232662"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2232654..2232662",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2244028",
                        "INSDInterval_to": "2244073"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2244197",
                        "INSDInterval_to": "2244503"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2244028..2244073,2244197..2244503)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.18.35"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2244504",
                    "INSDInterval_to": "2244510"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2244504..2244510",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2244534",
                    "INSDInterval_to": "2244542"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2244534..2244542",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2245935",
                    "INSDInterval_iscomp": {
                        "@value": "true"
                    },
                    "INSDInterval_to": "2245559"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "complement(2245559..2245935)",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Rpl26"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2263402",
                        "INSDInterval_to": "2263447"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2263532",
                        "INSDInterval_to": "2263834"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2263402..2263447,2263532..2263834)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.11.34"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2263836",
                    "INSDInterval_to": "2263842"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2263836..2263842",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2263866",
                    "INSDInterval_to": "2263874"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2263866..2263874",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2274582",
                        "INSDInterval_to": "2274627"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2274708",
                        "INSDInterval_to": "2275011"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2274582..2274627,2274708..2275011)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.10.33"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2275012",
                    "INSDInterval_to": "2275018"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2275012..2275018",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2275042",
                    "INSDInterval_to": "2275050"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2275042..2275050",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2276791",
                    "INSDInterval_to": "2277088"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2276791..2277088",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.6.32"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2289472",
                    "INSDInterval_to": "2289766"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2289472..2289766",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.17pg.31"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2298802",
                    "INSDInterval_to": "2299106"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2298802..2299106",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.5.30"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2300865",
                        "INSDInterval_to": "2300910"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2300995",
                        "INSDInterval_to": "2301298"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2300865..2300910,2300995..2301298)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.9.29"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2301299",
                    "INSDInterval_to": "2301305"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2301299..2301305",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2301329",
                    "INSDInterval_to": "2301337"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2301329..2301337",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2302969",
                    "INSDInterval_to": "2303285"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2302969..2303285",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.4.28"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2308440",
                        "INSDInterval_to": "2308485"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2308593",
                        "INSDInterval_to": "2308897"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2308440..2308485,2308593..2308897)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.16.27"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2308900",
                    "INSDInterval_to": "2308906"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2308900..2308906",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2308930",
                    "INSDInterval_to": "2308938"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2308930..2308938",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2325259",
                        "INSDInterval_to": "2325304"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2325424",
                        "INSDInterval_to": "2325722"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2325259..2325304,2325424..2325722)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.15pg.26"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2334582",
                        "INSDInterval_to": "2334627"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2334735",
                        "INSDInterval_to": "2335039"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2334582..2334627,2334735..2335039)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.14.25"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2335042",
                    "INSDInterval_to": "2335048"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2335042..2335048",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2335072",
                    "INSDInterval_to": "2335080"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2335072..2335080",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2346101",
                    "INSDInterval_to": "2346396"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2346101..2346396",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.13pg.24"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2351937",
                    "INSDInterval_to": "2352251"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2351937..2352251",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.3.23"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2354049",
                        "INSDInterval_to": "2354094"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2354179",
                        "INSDInterval_to": "2354482"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2354049..2354094,2354179..2354482)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.8.22"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2354483",
                    "INSDInterval_to": "2354489"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2354483..2354489",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2354513",
                    "INSDInterval_to": "2354521"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2354513..2354521",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2356255",
                    "INSDInterval_to": "2356544"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2356255..2356544",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.2.21"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2368573",
                        "INSDInterval_to": "2368618"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2368721",
                        "INSDInterval_to": "2369025"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2368573..2368618,2368721..2369025)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.12.20"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2369028",
                    "INSDInterval_to": "2369034"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2369028..2369034",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2369058",
                    "INSDInterval_to": "2369066"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2369058..2369066",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2383581",
                        "INSDInterval_to": "2383626"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2383752",
                        "INSDInterval_to": "2384029"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2383581..2383626,2383752..2384029)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.11pg.19"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2385341",
                        "INSDInterval_to": "2385386"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2385471",
                        "INSDInterval_to": "2385774"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2385341..2385386,2385471..2385774)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.7.18"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2385775",
                    "INSDInterval_to": "2385781"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2385775..2385781",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2385805",
                    "INSDInterval_to": "2385813"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2385805..2385813",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2388492",
                    "INSDInterval_to": "2388792"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2388492..2388792",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.6pg.17"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2401109",
                        "INSDInterval_to": "2401154"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2401277",
                        "INSDInterval_to": "2401575"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2401109..2401154,2401277..2401575)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.10pg.16"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2409028",
                        "INSDInterval_to": "2409073"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2409181",
                        "INSDInterval_to": "2409485"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2409028..2409073,2409181..2409485)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.9.15"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2409488",
                    "INSDInterval_to": "2409494"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2409488..2409494",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2409518",
                    "INSDInterval_to": "2409526"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2409518..2409526",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2415832",
                        "INSDInterval_to": "2415877"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2416002",
                        "INSDInterval_to": "2416289"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2415832..2415877,2416002..2416289)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.8pg.14"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2417532",
                        "INSDInterval_to": "2417577"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2417662",
                        "INSDInterval_to": "2417965"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2417532..2417577,2417662..2417965)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.5.13"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2417966",
                    "INSDInterval_to": "2417972"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2417966..2417972",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2417996",
                    "INSDInterval_to": "2418004"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2417996..2418004",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2423490",
                    "INSDInterval_to": "2423778"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2423490..2423778",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.4pg.12"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2436096",
                        "INSDInterval_to": "2436141"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2436264",
                        "INSDInterval_to": "2436563"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2436096..2436141,2436264..2436563)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "PG.1.11"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2440199",
                    "INSDInterval_to": "2441262"
                }
            }
        ],
        "INSDFeature_key": "gene",
        "INSDFeature_location": "2440199..2441262",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "gene",
                        "INSDQualifier_value": "Ankrd12"
                    },
                    {
                        "INSDQualifier_name": "pseudo"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2445300",
                        "INSDInterval_to": "2445345"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2445444",
                        "INSDInterval_to": "2445748"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2445300..2445345,2445444..2445748)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.7.10"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2445751",
                    "INSDInterval_to": "2445757"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2445751..2445757",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2445781",
                    "INSDInterval_to": "2445789"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2445781..2445789",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2457814",
                        "INSDInterval_to": "2457859"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2457969",
                        "INSDInterval_to": "2458262"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2457814..2457859,2457969..2458262)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.6pg.9"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2459639",
                        "INSDInterval_to": "2459684"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2459769",
                        "INSDInterval_to": "2460072"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2459639..2459684,2459769..2460072)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.3.8"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2460073",
                    "INSDInterval_to": "2460079"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2460073..2460079",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2460103",
                    "INSDInterval_to": "2460111"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2460103..2460111",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2467500",
                    "INSDInterval_to": "2467782"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2467500..2467782",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.5pg.7"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2473351",
                        "INSDInterval_to": "2473396"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2473504",
                        "INSDInterval_to": "2473808"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2473351..2473396,2473504..2473808)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.4.6"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2473811",
                    "INSDInterval_to": "2473817"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2473811..2473817",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2473841",
                    "INSDInterval_to": "2473849"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2473841..2473849",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2481014",
                    "INSDInterval_to": "2481291"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2481014..2481291",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.3pg.5"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2482556",
                        "INSDInterval_to": "2482601"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2482686",
                        "INSDInterval_to": "2482987"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2482556..2482601,2482686..2482987)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.2.4"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2482990",
                    "INSDInterval_to": "2482996"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2482990..2482996",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2483019",
                    "INSDInterval_to": "2483027"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2483019..2483027",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2492266",
                        "INSDInterval_to": "2492311"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2492446",
                        "INSDInterval_to": "2492750"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2492266..2492311,2492446..2492750)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.2.3"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2492753",
                    "INSDInterval_to": "2492759"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2492753..2492759",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2492783",
                    "INSDInterval_to": "2492791"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2492783..2492791",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2496705",
                    "INSDInterval_to": "2497015"
                }
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "2496705..2497015",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "Q52.1pg.2"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": [
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2497854",
                        "INSDInterval_to": "2497899"
                    },
                    {
                        "INSDInterval_accession": "BN000872.1",
                        "INSDInterval_from": "2498023",
                        "INSDInterval_to": "2498327"
                    }
                ]
            }
        ],
        "INSDFeature_key": "V_segment",
        "INSDFeature_location": "join(2497854..2497899,2498023..2498327)",
        "INSDFeature_operator": "join",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "7183.1pg.1"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2498330",
                    "INSDInterval_to": "2498336"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2498330..2498336",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS heptamer"
                    }
                ]
            }
        ]
    },
    {
        "INSDFeature_intervals": [
            {
                "INSDInterval": {
                    "INSDInterval_accession": "BN000872.1",
                    "INSDInterval_from": "2498360",
                    "INSDInterval_to": "2498368"
                }
            }
        ],
        "INSDFeature_key": "regulatory",
        "INSDFeature_location": "2498360..2498368",
        "INSDFeature_quals": [
            {
                "INSDQualifier": [
                    {
                        "INSDQualifier_name": "regulatory_class",
                        "INSDQualifier_value": "other"
                    },
                    {
                        "INSDQualifier_name": "note",
                        "INSDQualifier_value": "RSS nonamer"
                    }
                ]
            }
        ]
    }
]

buffer = {}
document = None
with io.open('/Users/lg/code/somatic/bootstrap/mouse_ighv_c57bl6_bn000872.json', 'r') as f:
    document = json.loads(f.read())
    for node in document:
        buffer[node['gene']] = node

node = None
for feature in INSDFeature:
    if 'INSDFeature_key' in feature:
        interval = feature['INSDFeature_intervals'][0]['INSDInterval']
        if isinstance(interval, list):
            interval = interval[-1]

        # print(to_json(feature))
        if feature['INSDFeature_key'] == 'V_segment':
            node = None
            name = feature['INSDFeature_quals'][0]['INSDQualifier'][-1]['INSDQualifier_value']
            if name in buffer:
                node = buffer[name]
                node['BN000872'] = {
                    'start': int(interval['INSDInterval_from']) - 1,
                    'end': int(interval['INSDInterval_to']),
                }
                node['BN000872']['length'] = node['BN000872']['end'] - node['BN000872']['start']

        elif feature['INSDFeature_key'] == 'regulatory':
            name = feature['INSDFeature_quals'][0]['INSDQualifier'][-1]['INSDQualifier_value']
            if 'heptamer' in name:
                name = 'heptamer'
            elif 'nonamer' in name:
                name = 'nonamer'
            if node is not None:
                node['BN000872'][name] = {
                    'start': int(interval['INSDInterval_from']) - 1,
                    'end': int(interval['INSDInterval_to'])
                }
  
            node['BN000872'][name]['length'] = node['BN000872'][name]['end'] - node['BN000872'][name]['start']

print(to_json(document))
