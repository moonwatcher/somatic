#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Somatic V/D/J recombination analyzer
# Author: Lior Galanti < lior.galanti@nyu.edu >
# NYU Center for Genetics and System Biology 2015
# { V } --> VDP5 --> VDN --> VDP3 --> { D } --> DJP5 --> DJN --> DJP3 --> { J }
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
import os
import logging
import re
import json
import uuid
import io
import hashlib
import math
import pymongo
import pickle
import numpy

from io import StringIO, BytesIO
from datetime import timedelta, datetime
from argparse import ArgumentParser
from subprocess import Popen, PIPE
from copy import deepcopy

from bson.objectid import ObjectId
from bson.binary import Binary
from pymongo.errors import BulkWriteError
from pymongo import MongoClient, DESCENDING, ASCENDING

log_levels = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

def verify_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)

def parse_match(match):
    if match is not None:
        return dict(((k.replace('_', ' '),v) for k,v in match.groupdict().items() if k and v))
    else:
        return None

def to_fastq(id, sequence, quality):
    if id:
        if sequence:
            if quality:
                return '@{}\n{}\n+\n{}'.format(id, sequence, quality)
            else:
                raise ValueError('FASTQ quality can not be empty')
        else:
            raise ValueError('FASTQ sequence can not be empty')
    else:
        raise ValueError('FASTQ id can not be empty')

def to_fasta(id, sequence, description, limit):
    if id and sequence:
        buffer = []
        if description:
            buffer.append('>{} {}'.format(id, description))
        else:
            buffer.append('>{}'.format(id))

        if limit:
            begin = 0
            end = 0
            length = len(sequence)
            while end < length:
                end = min(begin + limit, length)
                buffer.append(sequence[begin:end])
                begin = end
        else:
            buffer.append(sequence)
        return '\n'.join(buffer)
    else:
        raise ValueError('FASTA sequence must have an id and sequence')

def merge_base_quality(bi, bj, qi, qj):
    pi = phred_to_probability(qi)
    pj = phred_to_probability(qj)
    if bi == bj:
        pij = (pi * pj / 3.0) / (1.0 - pi - pj + ((4.0 * pi * pj) / 3.0))
    else:
        if qi > qj:
            pij = (pi * (1.0 - (pj / 3.0))) / (pi + pj - ((4.0 * pi * pj) / 3.0))
        else:
            pij = (pj * (1.0 - (pi / 3.0))) / (pj + pi - ((4.0 * pj * pi) / 3.0))

    qij = probability_to_phred(pij)
    if qij > 'I': qij = 'I'
    elif qij < '#': qij = '#'
    return qij

def phred_to_quality(code):
    return ord(code) - 33

def quality_to_phred(value):
    return chr(round(value) + 33)

def phred_to_probability(score):
    return pow(10, -(ord(score) - 33) / 10)

def probability_to_phred(probability):
    return chr(round(-10 * math.log(probability, 10)) + 33)

def expected_error(quality):
    return sum([ phred_to_probability(q) for q in quality ])

def to_csv(row):
    return ','.join([str(column) for column in row])

def to_document(node):
    if isinstance(node, set):
        return [ to_document(v) for v in node ]
    if isinstance(node, list):
        return [ to_document(v) for v in node ]
    elif isinstance(node, dict):
        return dict([ (k,to_document(v)) for k,v in node.items() ])
    elif isinstance(node, numpy.ndarray):
        if node.dtype == numpy.float:
            return [ float(v) for v in node ]
        elif node.dtype == numpy.int:
            return [ int(v) for v in node ]
        else:
            return node
    elif isinstance(node, Sample):
        return node.document
    elif isinstance(node, Artifact):
        return node.document
    elif isinstance(node, Gene):
        return node.document
    elif isinstance(node, Accession):
        return node.document
    elif isinstance(node, Sequence):
        return node.document
    elif isinstance(node, Study):
        return node.document
    elif isinstance(node, Pivot):
        return node.document
    elif isinstance(node, Request):
        return node.document
    elif isinstance(node, Hit):
        return node.document
    elif isinstance(node, Configuration):
        return node.document
    elif isinstance(node, Reference):
        return node.document
    else:
        return node

def to_json(node):
    def handler(o):
        result = None
        if isinstance(o, datetime):
            result = o.isoformat()
        if isinstance(o, ObjectId):
            result = str(o)
        return result
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4, default=handler)

def to_json_document(node):
    return to_json(to_document(node))

def simplify(value):
    if value: return value.lower()
    else: return None

def complement(nucleotide):
    complement = {
        'A': 'T',
        'B': 'V',
        'C': 'G',
        'D': 'H',
        'G': 'C',
        'H': 'D',
        'K': 'M',
        'M': 'K',
        'N': 'N',
        'R': 'Y',
        'S': 'S',
        'T': 'A',
        'V': 'B',
        'W': 'W',
        'Y': 'R',
        'a': 't',
        'b': 'v',
        'c': 'g',
        'd': 'h',
        'g': 'c',
        'h': 'd',
        'k': 'm',
        'm': 'k',
        'n': 'n',
        'r': 'y',
        's': 's',
        't': 'a',
        'v': 'b',
        'w': 'w',
        'y': 'r',
    }
    return ''.join([ complement[n] for n in nucleotide ][::-1])

def merge(this, other):
    result = other
    if this is not None:
        result = this
        if other is not None:
            if isinstance(this, dict):
                if isinstance(other, dict):
                    for k,v in other.items():
                        if k in result:
                            result[k] = merge(result[k], v)
                        else:
                            result[k] = v
                else:
                    raise ValueError('incompatible structure')

            elif isinstance(this, list):
                if isinstance(other, list):
                    result.extend(other)
                else:
                    raise ValueError('incompatible structure')
            else:
                result = other
    return result

def merge_configuration_node(node, other):
    result = other
    if node is not None:
        result = node
        if other is not None:
            if isinstance(node, dict):
                if isinstance(other, dict):
                    for k,v in other.items():
                        if k in result:
                            result[k] = merge_configuration_node(result[k], v)
                        else:
                            result[k] = v
                else:
                    raise ValueError('invalid configuration structure')

            elif isinstance(node, list):
                if isinstance(other, list):
                    result.extend(other)
                else:
                    raise ValueError('invalid configuration structure')
            else:
                result = other
    return result

def check_configuration_node(element):
    result = None
    if isinstance(element, dict):
        if 'enabled' not in element or element['enabled']:
            if 'enabled' in element:
                del element['enabled']

            for k in list(element.keys()):
                element[k] = check_configuration_node(element[k])
            result = element

    elif isinstance(element, list):
        result = []
        for o in element:
            checked = check_configuration_node(o)
            if checked is not None:
                result.append(checked)
    else:
        result = element

    return result

def hamming(one, two):
    distance = 0
    for o, t in zip(one, two):
      if o != t:
          distance += 1
    return distance

def sequence_diff(one, two):
    diff = []
    for o, t in zip(one, two):
        if o == t: diff.append('-')
        else: diff.append('*')
    return ''.join(diff)

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

class Command(object):
    def __init__(self, pipeline, name, node=None):
        self.log = logging.getLogger('Command')
        self.pipeline = pipeline
        self.node = {
            'stdin': PIPE,
            'stdout': PIPE,
            'stderr': PIPE,
            'cwd': None,
            'environment': None,
            'binary': None
        }
        self.merge(deepcopy(self.configuration.command[name]))
        self.merge(node)

    def merge(self, node):
        self.node = merge(self.node, node)

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def binary(self):
        return self.node['binary']

    @property
    def parameter(self):
        if 'parameter' not in self.node:
            self.node['parameter'] = {}
        return self.node['parameter']

    @property
    def environment(self):
        return None if 'environment' not in self.node else self.node['environment']

    @property
    def cwd(self):
        return None if 'cwd' not in self.node else self.node['cwd']

    @cwd.setter
    def cwd(self, value):
        self.node['cwd'] = value

    @property
    def stdin(self):
        return self.node['stdin']

    @stdin.setter
    def stdin(self, value):
        self.node['stdin'] = value

    @property
    def stdout(self):
        return self.node['stdout']

    @stdout.setter
    def stdout(self, value):
        self.node['stdout'] = value

    @property
    def stderr(self):
        return self.node['stderr']

    @stderr.setter
    def stderr(self, value):
        self.node['stderr'] = value

    @property
    def args(self):
        if 'args' not in self.node:
            self.node['args'] = [ self.binary ]
            if 'sub command' in self.node:
                self.node['args'].append(self.node['sub command'])

            if self.parameter:
                for k,v in self.parameter.items():
                    self.node['args'].append(k)
                    if v is not None:
                        self.node['args'].append(v)
        return self.node['args']

    @property
    def process(self):
        if 'process' not in self.node:
            self.node['process'] = Popen(
                args=self.args,
                env=self.environment,
                cwd=self.cwd,
                stdin=self.stdin,
                stdout=self.stdout,
                stderr=self.stderr,
            )
        return self.node['process']

class Configuration(object):
    def __init__(self, node=None):
        self.log = logging.getLogger('Configuration')
        self.node = {
            'command': None,
            'interface': None,
            'constant': None,
            'template': None,
            'preset': None,
            'profile': None,
            'expression': {},
            'diagram': {},
            'enumeration': None,
            'rss': None,
            'database': None,
            'extension': None,
            'reference': None,
        }
        self.selected = {
            'organism': None,
            'strain': None,
            'chain': None,
            'region': None,
        }
        self.merge(node)

        path = os.path.expanduser(os.path.expandvars('~/.somatic/core.json'))
        with io.open(path, 'rb') as file:
            self.merge(json.loads(file.read().decode('utf8')))

        path = os.path.expanduser(os.path.expandvars('~/.somatic/setting.json'))
        with io.open(path, 'rb') as file:
            self.merge(json.loads(file.read().decode('utf8')))

        self.load_cache()
        self.load_reference()
        self.load_organism()
        self.load_expression()
        self.load_preset()
        self.load_diagram()
        self.load_interface()
        self.load_stop_codon()
        self.load_reverse_complement()

    def load_cache(self):
        self.cache['path'] = os.path.expanduser(os.path.expandvars(self.cache['path']))
        self.cache['resource directory'] = os.path.join(self.cache['path'], 'resource')
        self.cache['reference directory'] = os.path.join(self.cache['path'], 'reference')
        self.cache['download directory'] = os.path.join(self.cache['path'], 'download')

        verify_directory(self.cache['path'])
        verify_directory(self.cache['resource directory'])
        verify_directory(self.cache['download directory'])
        verify_directory(self.cache['reference directory'])

    def load_reference(self):
        for key, reference in self.reference.items():
            reference['key'] = key
            basename = os.path.split(reference['remote'])[1]
            filename, extension = os.path.splitext(basename)
            extension = extension.strip('.')
            if extension in ['gz', 'bz2']:
                reference['path'] = os.path.join(self.cache['reference directory'], '{}.pickled'.format(filename))
                reference['local'] = os.path.join(self.cache['download directory'], filename)
            else:
                reference['path'] = os.path.join(self.cache['reference directory'], '{}.pickled'.format(basename))
                reference['local'] = os.path.join(self.cache['download directory'], basename)

    def load_organism(self):
        for o, organism in self.organism.items():
            organism['key'] = o
            organism['resource directory'] = os.path.join(self.cache['resource directory'], organism['key'].replace(' ', '_'))
            if 'strain' in organism:
                for s, strain in organism['strain'].items():
                    if strain is not None:
                        strain['key'] = s
                        strain['resource directory'] = os.path.join(organism['resource directory'], strain['key'].replace(' ', '_'))
                    if 'chain' in strain:
                        for c, chain in strain['chain'].items():
                            if chain is not None:
                                chain['key'] = c
                                chain['resource directory'] = os.path.join(strain['resource directory'], chain['key'].replace(' ', '_'))
                                chain['bwa directory'] = os.path.join(chain['resource directory'], 'bwa')
                                chain['igblast directory'] = os.path.join(chain['resource directory'], 'igblast')
                                chain['igblast database directory'] = os.path.join(chain['igblast directory'], 'database')
                                if 'region' in chain:
                                    for r, region in chain['region'].items():
                                        if region is not None:
                                            region['key'] = r
                                            if region['type'] == 'gene segment':
                                                region['igblast fasta path'] = os.path.join(chain['igblast database directory'], '{}.fa'.format(region['key']))
                                                region['bwa fasta path'] = os.path.join(chain['bwa directory'], '{}.fa'.format(region['key']))

    def load_diagram(self):
        self.diagram['prototype'] = {
            'name': {
                'format': lambda x: x,
                'title': 'name',
                'width': 'auto',
                'value': 'gene'
            },
            'strain': {
                'format': lambda x: x,
                'title': 'strain',
                'width': 'auto',
                'value': 'strain',
            },
            'phred': {
                'format': lambda x: '' if x is None else '{:.3}'.format(x),
                'title': 'Q',
                'width': 'auto',
                'value': 'average phread',
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

    def load_preset(self):
        for key in list(self.template['study']['feature'].keys()):
            feature = deepcopy(self.template['study']['default']['feature'])
            for column in list(feature['column'].keys()):
                feature['column'][column] = merge_configuration_node(deepcopy(self.template['study']['default']['column']), feature['column'][column])
                feature['column'][column]['key'] = column
            feature['key'] = key
            feature = merge_configuration_node(feature, self.template['study']['feature'][key])
            self.template['study']['feature'][key] = feature

        for name in list(self.preset.keys()):
            if 'study' in self.preset[name]:
                preset = deepcopy(self.template['study']['default']['preset'])
                preset['key'] = name
                preset = merge_configuration_node(preset, self.preset[name]['study'])

                for key in list(preset['feature'].keys()):
                    feature = deepcopy(self.template['study']['feature'][key])
                    feature = merge_configuration_node(feature, preset['feature'][key])
                    preset['feature'][key] = feature

                row = []
                for feature_key in preset['order']:
                    if feature_key in preset['feature']:
                        feature = { 'key': feature_key, 'column': [] }
                        for column_key in preset['feature'][feature_key]['order']:
                            if column_key in preset['feature'][feature_key]['column']:
                                feature['column'].append(preset['feature'][feature_key]['column'][column_key])
                        row.append(feature)
                preset['row'] = check_configuration_node(row)
                del preset['feature']
                del preset['order']
                self.preset[name]['study'] = preset

    def load_expression(self):
        self.expression['cigar'] = re.compile('\*|([0-9]+[MIDNSHPX=])+')
        self.expression['ncbi accession url'] = 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&dopt=gbc_xml&val={}'
        self.expression['igblast reversed query'] = '# Note that your query represents the minus strand of a V gene and has been converted to the plus strand.'
        self.expression['nucleotide sequence'] = re.compile('^[ACGTRYKMSWBdhVN]+$', re.IGNORECASE)
        self.expression['sam record'] = re.compile(
            r"""^
            (?P<QNAME>[!-?A-~]{1,254})\t
            (?P<FLAG>[0-9]+)\t
            (?P<RNAME>\*|[!-()+-<>-~][!-~]*)\t
            (?P<POS>[0-9]+)\t
            (?P<MAPQ>[0-9]+)\t
            (?P<CIGAR>\*|([0-9]+[MIDNSHPX=])+)\t
            (?P<RNEXT>\*|=|[!-()+-<>-~][!-~]*)\t
            (?P<PNEXT>[0-9]+)\t
            (?P<TLEN>[-+]?[0-9]+)\t
            (?P<SEQ>\*|[A-Za-z=.]+)\t
            (?P<QUAL>[!-~]+)
            (?P<OPT>(\t([A-Za-z][A-Za-z0-9]:[AifZHB]:[^:]+))+)?
            $""",
            re.VERBOSE
        )
        self.expression['sam optional field'] = re.compile('(?P<TAG>[A-Za-z][A-Za-z0-9]):(?P<TYPE>[AifZHB]):(?P<VALUE>[^:]+)')
        self.expression['sam optional field value'] = {
            'A': re.compile('^[!-~]$'),
            'i': re.compile('^[-+]?[0-9]+$'),
            'f': re.compile('^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'),
            'Z': re.compile('^[ !-~]+$'),
            'H': re.compile('^[0-9A-F]+$'),
            'B': re.compile('^[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+$'),
        }
        self.expression['blat hit'] = re.compile(
            r"""
            (?P<match>[^\t]+)\t
            (?P<mismatch>[^\t]+)\t
            (?P<repeat_match>[^\t]+)\t
            (?P<n_count>[^\t]+)\t
            (?P<inserts_in_query>[^\t]+)\t
            (?P<inserted_base_in_query>[^\t]+)\t
            (?P<inserts_in_target>[^\t]+)\t
            (?P<inserted_base_in_target>[^\t]+)\t
            (?P<query_strand>[^\t]+)\t
            (?P<query_name>[^\t]+)\t
            (?P<query_size>[^\t]+)\t
            (?P<query_start>[^\t]+)\t
            (?P<query_end>[^\t]+)\t
            (?P<target_name>[^\t]+)\t
            (?P<target_size>[^\t]+)\t
            (?P<target_start>[^\t]+)\t
            (?P<target_end>[^\t]+)\t
            (?P<block_count>[^\t]+)\t
            (?P<block_size>[^\t]+)\t
            (?P<query_block_start>[^\t]+)\t
            (?P<target_block_start>[^\t]+)
            """,
            re.VERBOSE
        )
        self.expression['igblast hit'] = re.compile(
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
        )

    def load_interface(self):
        self.interface['prototype']['preset']['parameter']['choices'] = list(self.preset.keys())
        self.interface['prototype']['preset']['parameter']['help'] = 'one of {}'.format(', '.join(sorted(self.preset.keys())))

        self.interface['prototype']['profile']['parameter']['choices'] = list(self.profile.keys())
        self.interface['prototype']['profile']['parameter']['help'] = 'one of: {}'.format(', '.join(sorted(self.profile.keys())))

        self.interface['implementation'] = {}
        for action in self.interface['section']['action']:
            self.interface['implementation'][action['instruction']['name']] = action['implementation']

    def load_reverse_complement(self):
        self.enumeration['reverse complement'] = {}
        for k,v in self.enumeration['iupac nucleic acid notation'].items():
            self.enumeration['reverse complement'][k] = v['reverse']

    def load_stop_codon(self):
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

        seed = [ k for k,v in self.enumeration['nucleic to amino'].items() if v == '*' ]
        feature = {
            'codon': dict([ (c, { 'motif': [] }) for c in seed ]),
            'space': set()
        }

        for triplet, codon in feature['codon'].items():
            for nucleotide in triplet:
                motif = []
                for k,n in self.enumeration['iupac nucleic acid notation'].items():
                    if nucleotide in n['option']:
                        motif.append(k)
                codon['motif'].append(motif)

            codon['possible'] = expand(codon['motif'])
            feature['space'] |= set(codon['possible'])
        self.enumeration['stop codon repertoire'] = feature['space']

    def merge(self, other):
        self.node = merge_configuration_node(self.node, other)

    def select(self, instruction):
        if 'organism' in instruction:
            self.selected['organism'] = self.organism[instruction['organism']]
            if 'strain' in instruction:
                self.selected['strain'] = self.selected['organism']['strain'][instruction['strain']]
                if 'chain' in instruction:
                    self.selected['chain'] = self.selected['strain']['chain'][instruction['chain']]
                    self.selected['region'] = self.selected['chain']['region']

    @property
    def document(self):
        return to_document(self.node)

    @property
    def json(self):
        return to_json_document(self.node)

    @property
    def command(self):
        return self.node['command']

    @property
    def interface(self):
        return self.node['interface']

    @property
    def constant(self):
        return self.node['constant']

    @property
    def template(self):
        return self.node['template']

    @property
    def preset(self):
        return self.node['preset']

    @property
    def profile(self):
        return self.node['profile']

    @property
    def expression(self):
        return self.node['expression']

    @property
    def diagram(self):
        return self.node['diagram']

    @property
    def enumeration(self):
        return self.node['enumeration']

    @property
    def reference(self):
        return self.node['reference']

    @property
    def rss(self):
        return self.node['rss']

    @property
    def organism(self):
        return self.node['organism']

    @property
    def database(self):
        return self.node['database'][self.system['host']]

    @property
    def system(self):
        return self.node['system']

    @property
    def extension(self):
        return self.node['extension']

    @property
    def table(self):
        return self.node['table']

    @property
    def cache(self):
        return self.node['cache']

class Blat(object):
    def __init__(self, pipeline, instruction):
        self.log = logging.getLogger('Blat')
        self.pipeline = pipeline
        self.instruction = instruction
        self.reference = self.pipeline.resolver.reference_fetch(self.configuration.selected['chain']['reference'])
        self.fasta = []
        for record in self.instruction['record'].values():
            record['flanking'] = record['gene'].flanking_accession_query(instruction['flank'])
            self.fasta.append(
                to_fasta(
                    record['flanking']['id'],
                    record['flanking']['sequence'].nucleotide,
                    None,
                    self.configuration.constant['fasta line length']
                )
            )
        self.fasta = '\n'.join(self.fasta)

    @property
    def configuration(self):
        return self.pipeline.configuration

    def orientation(self, query, hit):
        # if the query is on the opposite strand from the target, reverse the query
        orientation = { 'start': None, 'end': None, 'flanking': None }
        if hit['query strand']:
            orientation['flanking'] = query['flanking']['sequence']
            orientation['end'] = query['flanking']['end']
            orientation['start'] = query['flanking']['start']
        else:
            orientation['flanking'] = query['flanking']['sequence'].reversed
            orientation['end'] = query['flanking']['flank length'] - query['flanking']['start']
            orientation['start'] = query['flanking']['flank length'] - query['flanking']['end']
        return orientation

    def search(self):
        strain = self.configuration.selected['strain']
        chain = self.configuration.selected['chain']
        command = Command(self.pipeline, 'blat', chain['blat'])
        command.args.insert(1, self.reference.local)
        command.args.insert(2, 'stdin')
        command.args.insert(3, 'stdout')
        output, error = command.process.communicate(input=self.fasta.encode('utf8'))
        if output:
            buffer = StringIO(output.decode('utf8'))
            for line in buffer:
                hit = parse_match(self.configuration.expression['blat hit'].search(line.strip()))
                if hit:
                    query = self.instruction['record'][hit['query name']]
                    self.parse_flanked_hit(query, hit)

            # for every query sort the results by the objective score and than the flanked score
            for query in self.instruction['record'].values():
                if 'hit' in query:
                    self.crop(query)
                    query['hit'].sort(key=lambda x: x['flanked']['identical'], reverse=True)
                    query['hit'].sort(key=lambda x: x['objective']['identical'], reverse=True)
                    query['hit'].sort(key=lambda x: x['flanked']['score'], reverse=True)
                    query['hit'].sort(key=lambda x: x['objective']['score'], reverse=True)
                    self.summarize(query)

    def parse_flanked_hit(self, query, flanked):
        flanked['query strand'] = flanked['query strand'] == '+'
        
        for k in [
            'block size',
            'query block start',
            'target block start',
        ]:
            if k in flanked:
                flanked[k] = [ int(v) for v in flanked[k].split(',') if v ]
            
        for k in [
            'block count',
            'inserted base in query',
            'inserted base in target',
            'inserts in query',
            'inserts in target',
            'match',
            'mismatch',
            'n count',
            'query end',
            'query size',
            'query start',
            'repeat match',
            'target end',
            'target size',
            'target start',
        ]:
            if k in flanked:
                flanked[k] = int(flanked[k])
            
        flanked['block'] = []
        orientation = self.orientation(query, flanked)
        record = self.reference.record[flanked['target name']]
        flanked['reference file sha1'] = self.reference.sha1
        flanked['reference record sha1'] = record['sha1']
        for i in range(flanked['block count']):
            block = {
                'size': flanked['block size'][i],
                'query start': flanked['query block start'][i],
                'target start': flanked['target block start'][i],
            }
            block['target sequence'] = record['sequence'].nucleotide[block['target start']:block['target start'] + block['size']]
            block['query sequence'] = orientation['flanking'].nucleotide[block['query start']:block['query start'] + block['size']]
            flanked['block'].append(block)
        del flanked['block size']
        del flanked['query block start']
        del flanked['target block start']
        del flanked['query name']

        if 'hit' not in query: query['hit'] = []
        query['hit'].append({'flanked': self.normalize_hit(flanked, orientation)})

    def normalize_hit(self, hit, orientation):
        # construct a consensus sequence for both query and target
        position = { 'query': hit['query start'], 'target': hit['target start'] }
        hit['query sequence'] = ''
        hit['target sequence'] = ''

        for block in hit['block']:
            if block['query start'] > position['query']:
                gap = block['query start'] - position['query']
                hit['query sequence'] += orientation['flanking'].nucleotide[position['query']:block['query start']]
                position['query'] += gap
                hit['target sequence'] += '-' * gap
                
            if block['target start'] > position['target']:
                record = self.reference.record_by_sha1(hit['reference record sha1'])
                gap = block['target start'] - position['target']
                hit['target sequence'] += record['sequence'].nucleotide[position['target']:position['target'] + gap]
                position['target'] += gap
                hit['query sequence'] += '-' * gap
                
            hit['query sequence'] += block['query sequence']
            position['query'] += block['size']
            
            hit['target sequence'] += block['target sequence']
            position['target'] += block['size']
            
        # check the result
        hit['check'] = []
        for i in range(len(hit['query sequence'])):
            q = hit['query sequence'][i]
            t = hit['target sequence'][i]
            if q == '-' or t == '-' or q == t:
                hit['check'].append('-')
            else:
                hit['check'].append('*')
        hit['check'] = ''.join(hit['check'])

        hit['query length'] = hit['query end'] - hit['query start']
        hit['target length'] = hit['target end'] - hit['target start']
        hit['score'] = hit['match'] + (hit['repeat match'] / 2) - hit['mismatch'] - hit['inserts in query'] - hit['inserts in target']
        
        millibad = 0
        if hit['query length'] > 0 and hit['target length'] > 0:
            diff = max(hit['query length'] - hit['target length'], 0)
            total = hit['match'] + hit['repeat match'] + hit['mismatch']
            if total != 0:
                millibad = (1000 * (hit['mismatch'] + hit['inserts in query'] + round(3 * math.log(diff + 1)))) / total
        hit['identical'] = 100.0 - millibad * 0.1
        return hit

    def crop(self, query):
        if 'hit' in query:
            for hit in query['hit']:
                orientation = self.orientation(query, hit['flanked'])
                record = self.reference.record_by_sha1(hit['flanked']['reference record sha1'])

                # extract both query and target sequence for the contiguous blocks
                for block in hit['flanked']['block']:
                    block['target sequence'] = record['sequence'].nucleotide[block['target start']:block['target start'] + block['size']]
                    block['query sequence'] = orientation['flanking'].nucleotide[block['query start']:block['query start'] + block['size']]
                    
                hit['objective'] = {
                    'block': [],
                    'score': 0,
                    'identical': 0.0,
                    'match': 0,
                    'mismatch': 0,
                    'block count': 0,
                    'reference file sha1': hit['flanked']['reference file sha1'],
                    'reference record sha1': hit['flanked']['reference record sha1'],
                    'query start': orientation['start'],
                    'query end': orientation['end'],
                    'query strand': hit['flanked']['query strand'],
                    'repeat match': hit['flanked']['repeat match'],
                    'inserted base in query': 0,
                    'inserted base in target': 0,
                    'inserts in query': 0,
                    'inserts in target': 0,
                }

                # project the contiguous blocks on the objective region
                for block in hit['flanked']['block']:
                    o = {
                        'match': 0,
                        'mismatch': 0,
                        'query start': max(block['query start'], orientation['start']),
                        'query end': min(block['query start'] + block['size'], orientation['end'])
                    }
                    if o['query end'] > o['query start']:
                        o['size'] = o['query end'] - o['query start']
                        o['query sequence'] = orientation['flanking'].nucleotide[o['query start']:o['query end']]
                        o['target start'] = o['query start'] + block['target start'] - block['query start']
                        o['target end'] = o['target start'] + o['size']
                        o['target sequence'] = record['sequence'].nucleotide[o['target start']:o['target start'] + o['size']]

                        for i in range(o['size']):
                            if o['target sequence'][i] == o['query sequence'][i]:
                                o['match'] += 1
                            else:
                                o['mismatch'] += 1
                        hit['objective']['match'] += o['match']
                        hit['objective']['mismatch'] += o['mismatch']
                        hit['objective']['block count'] += 1
                        hit['objective']['block'].append(o)

                if hit['objective']['block']:
                    hit['objective']['target start'] = min([ b['target start'] for b in hit['objective']['block'] ])
                    hit['objective']['target end'] = max([ b['target end'] for b in hit['objective']['block'] ])
                    
                    # infer the gaps
                    for i in range(len(hit['objective']['block']) - 1):
                        gap = hit['objective']['block'][i + 1]['query start'] - hit['objective']['block'][i]['query end']
                        if gap:
                            hit['objective']['inserts in query'] += 1
                            hit['objective']['inserted base in query'] += gap

                        gap = hit['objective']['block'][i + 1]['target start'] - hit['objective']['block'][i]['target end']
                        if gap:
                            hit['objective']['inserts in target'] += 1
                            hit['objective']['inserted base in target'] += gap
                    self.normalize_hit(hit['objective'], orientation)

    def summarize(self, query):
        if query is not None:
            query['summary'] = []
            for hit in query['hit']:
                if ('objective' in hit and len(hit['objective']['block']) > 0 and
                    hit['objective']['score'] == query['hit'][0]['objective']['score'] and
                    abs(hit['flanked']['score'] - query['hit'][0]['flanked']['score']) < 0.5 and
                    hit['objective']['identical'] == query['hit'][0]['objective']['identical'] and
                    abs(hit['flanked']['identical'] - query['hit'][0]['flanked']['identical']) < 0.5):

                    match = deepcopy(hit['objective'])
                    match['flanked score'] = hit['flanked']['score']
                    match['flanked identical'] = hit['flanked']['identical']
                    match['flanked query size'] = hit['flanked']['query size']
                    match['flanked target length'] = hit['flanked']['target length']
                    query['summary'].append(match)

            if not query['summary']:
                del query['summary']

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
    def valid(self):
        return self.length > 0

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def document(self):
        document = to_document(self.node)
        for key in [
            'codon',
            'phred',
        ]:
            if key in document:
                del document[key]
        return document

    def reset(self):
        if self.node is None:
            self.node = { 'read frame': 0, 'strand': True }
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

    def crop(self, start, end=None):
        sequence = None
        start = 0 if start is None else max(start, 0)
        end = self.length if end is None else min(end, self.length)
        if end > start:
            cropped = {
                'nucleotide': self.nucleotide[start:end],
                'strand': self.strand
            }
            if self.read_frame is not None:
                cropped['read frame'] = (3 - (start - self.read_frame)%3)%3

            if 'quality' in self.node:
                cropped['quality'] = self.quality[start:end]
            sequence = Sequence(self.pipeline, cropped)
        return sequence


        # sequence = None
        # if start is not None:
        #     start = max(start, 0)
        #     end = self.length if end is None else min(end, self.length)
        #     if end > start:
        #         cropped = {
        #             'nucleotide': self.nucleotide[start:end],
        #             'strand': self.strand
        #         }
        #         if 'read frame' in self.node:
        #             cropped['read frame'] = (3 - (start - self.read_frame)%3)%3

        #         if 'quality' in self.node:
        #             cropped['quality'] = self.quality[start:end]
        #         sequence = Sequence(self.pipeline, cropped)
        # return sequence

    @property
    def read_frame(self):
        if 'read frame' in self.node:
            return self.node['read frame']
        else:
            return None

    @read_frame.setter
    def read_frame(self, value):
        if value is None:
            if 'read frame' in self.node:
                del self.node['read frame']
        else:
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
            return ''

    @nucleotide.setter
    def nucleotide(self, value):
        if value is None:
            if 'nucleotide' in self.node:
                del self.node['nucleotide']
        else:
            self.node['nucleotide'] = value
        self.reset()

    @property
    def quality(self):
        if 'quality' in self.node:
            return self.node['quality']
        else:
            return ''

    @quality.setter
    def quality(self, value):
        if value is None:
            if 'quality' in self.node:
                del self.node['quality']
        else:
            self.node['quality'] = value
        self.reset()

    @property
    def codon(self):
        if 'codon' not in self.node and self.nucleotide:
            codon = self.codon_at_frame(self.read_frame)
            if codon:
                self.node['codon'] = codon
        return None if 'codon' not in self.node else self.node['codon']

    @property
    def classified(self):
        classified = None
        if self.codon:
            classified = []
            for codon in self.codon:
                acid = self.configuration.enumeration['iupac amino acid notation'][codon]
                if 'class code' in acid:
                    classified.append(acid['class code'])
                else:
                    classified.append('X')
            classified = ''.join(classified)
        return classified

    @property
    def phred(self):
        if 'phred' not in self.node and self.quality:
            self.node['phred'] = [ (ord(c) - 33) for c in self.quality ]
        return self.node['phred']

    def codon_at_frame(self, read_frame):
        start = read_frame
        end = start + 3
        codon = []
        while(not end > self.length):
            acid = self.nucleotide[start:end]
            try:
                codon.append(self.configuration.enumeration['nucleic to amino'][acid])
            except KeyError as e:
                self.log.debug('could not resolve %s triplet to an amino acid', e)
                if acid in self.configuration.enumeration['stop codon repertoire']:
                    self.log.debug('%s is potentially a stop codon', acid)
                    codon.append('*')
                else:
                    codon.append('X')
            start = end
            end = start + 3
        return ''.join(codon)

    @property
    def length(self):
        return len(self.nucleotide)

    @property
    def reversed(self):
        if self._reversed is None and self.nucleotide:
            reversed = {
                'nucleotide': ''.join([ self.configuration.enumeration['reverse complement'][b] for b in list(self.nucleotide) ][::-1]),
                'strand': not self.strand
            }
            if 'read frame' in self.node:
                reversed['read frame'] = (self.length - self.read_frame) % 3

            if 'quality' in self.node:
                reversed['quality'] = self.quality[::-1]
                
            self._reversed = Sequence(self.pipeline, reversed)
        return self._reversed

class Accession(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Accession')
        self.pipeline = pipeline
        self.node = node
        
        if self.node is None:
            self.node = {
                'head': {},
                'body': {},
            }
            
        if 'sequence' in self.body and isinstance(self.body['sequence'], dict):
            self.body['sequence'] = Sequence(self.pipeline, self.body['sequence'])
        else:
            self.body['sequence'] = Sequence(self.pipeline)

    def __str__(self):
        buffer = []
        buffer.append(self.id if self.id else 'unknown id')
        buffer.append(self.stain if self.strain else 'unknown strain')
        buffer.append('{} nucleotides'.format(self.length))
        return '[ {} ]'.format(', '.join(buffer))

    @property
    def valid(self):
        return self.id is not None and self.sequence.valid

    @property
    def document(self):
        return to_document(self.node)

    @property
    def fasta(self):
        return self.to_fasta(self.configuration.constant['fasta line length'])

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
        if 'id' in self.head:
            return self.head['id']
        else:
            return None

    @id.setter
    def id(self, value):
        self.head['id'] = value
        if self.head['id'] == '':
            del self.head['id']

    @property
    def sequence(self):
        return self.body['sequence']

    @property
    def length(self):
        return self.sequence.length

    @property
    def description(self):
        if 'description' in self.head:
            return self.head['description']
        else:
            return None

    @property
    def strain(self):
        if 'strain' in self.head:
            return self.head['strain']
        else:
            return None

    @strain.setter
    def strain(self, value):
        self.head['strain'] = value
        if self.head['strain'] == '':
            del self.head['strain']

    def to_fasta(self, limit):
        return to_fasta(self.id, self.sequence.nucleotide, self.description, limit)

    def to_json(self):
        print(to_json_document(self))

class Artifact(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Artifact')
        self.pipeline = pipeline
        self.node = node
        
        if self.node is None:
            self.node = {
                'head': { },
                'body': { 'accession strand': True }
            }
            
        if 'sequence' in self.body:
            if isinstance(self.body['sequence'], dict):
                self.body['sequence'] = Sequence(self.pipeline, self.body['sequence'])
        else:
            self.body['sequence'] = Sequence(self.pipeline)

    def __str__(self):
        buffer = []
        buffer.append(self.id if self.id else 'unknown id')
        buffer.append(self.region if self.region else 'unknown region')
        buffer.append(self.strain if self.strain else 'unknown strain')
        if 'accession' in self.head: buffer.append(self.head['accession'])
        if 'accession strand' in self.body: buffer.append('+' if self.accession_strand else '-')
        buffer.append('{}bp'.format(self.length))
        return '[ {} ]'.format(', '.join(buffer))

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def head(self):
        return self.node['head']

    @property
    def body(self):
        return self.node['body']

    @property
    def document(self):
        return to_document(self.node)

    @property
    def fasta(self):
        return to_fasta(self.id, self.sequence.nucleotide, None, self.configuration.constant['fasta line length'])

    @property
    def valid(self):
        return self.id is not None and self.sequence.valid

    @property
    def id(self):
        if 'id' in self.head:
            return self.head['id']
        else:
            return None

    @property
    def organism(self):
        if 'organism name' in self.head:
            return self.head['organism name']
        else:
            return None

    @property
    def strain(self):
        if 'strain' in self.head:
            return self.head['strain']
        else:
            return None

    @property
    def region(self):
        if 'region' in self.head:
            return self.head['region']
        else:
            return None

    @property
    def accession(self):
        if 'accession' in self.head:
            return self.head['accession']
        else:
            return None

    @property
    def aligned(self):
        if 'aligned' in self.head:
            return self.head['aligned']
        else:
            return False

    @property
    def alignment(self):
        if 'alignment' in self.body:
            return self.body['alignment']
        else:
            return None
    
    @property
    def reference_start(self):
        if 'reference start' in self.body:
            return self.body['reference start']
        else:
            return None
    
    @property
    def reference_end(self):
        if 'reference end' in self.body:
            return self.body['reference end']
        else:
            return None
    
    @property
    def reference_strand(self):
        if 'reference strand' in self.body:
            return self.body['reference strand']
        else:
            return None
    
    @property
    def accession_start(self):
        if 'accession start' in self.body:
            return self.body['accession start']
        else:
            return None

    @property
    def accession_end(self):
        if 'accession end' in self.body:
            return self.body['accession end']
        else:
            return None

    @property
    def accession_strand(self):
        return self.body['accession strand']

    @property
    def length(self):
        return self.sequence.length

    @property
    def sequence(self):
        return self.body['sequence']

    def flanking_accession_query(self, flank):
        query = None
        if self.accession:
            accession = self.pipeline.resolver.accession_fetch(self.accession)
            if accession is not None:
                query = {}
                query['id'] = self.id
                query['accession start'] = self.accession_start
                query['accession end'] = self.accession_end
                query['flank start'] = max(self.accession_start - flank, 0)
                query['flank end'] = min(self.accession_end + flank, accession.length)
                query['flank length'] = query['flank end'] - query['flank start']
                if self.accession_strand:
                    query['start'] = query['accession start'] - query['flank start']
                    query['end'] = query['accession end'] - query['flank start']
                    query['sequence'] = accession.sequence.crop(query['flank start'], query['flank end'])
                else:
                    query['start'] = query['flank end'] - query['accession end']
                    query['end'] = query['flank end'] - query['accession start']
                    query['sequence'] = accession.sequence.crop(query['flank start'], query['flank end']).reversed
                query['length'] = query['end'] - query['start']
        return query

    def flanking_reference_query(self, flank):
        query = None
        if self.aligned:
            reference = self.pipeline.resolver.reference_fetch(self.body['reference file sha1'])
            if reference is not None:
                record = reference.record_by_sha1(self.body['reference record sha1'])
                query = {}
                query['id'] = self.id
                query['reference start'] = self.body['reference start']
                query['reference end'] = self.body['reference end']
                query['reference strand'] = self.body['reference strand']
                query['flank start'] = max(query['reference start'] - flank, 0)
                query['flank end'] = min(query['reference end'] + flank, record['sequence'].length)
                query['flank length'] = query['flank end'] - query['flank start']
                if self.body['reference strand']:
                    query['start'] = query['reference start'] - query['flank start']
                    query['end'] = query['reference end'] - query['flank start']
                    query['sequence'] = record['sequence'].crop(query['flank start'], query['flank end'])
                else:
                    query['start'] = query['flank end'] - query['reference end']
                    query['end'] = query['flank end'] - query['reference start']
                    query['sequence'] = record['sequence'].crop(query['flank start'], query['flank end']).reversed
                query['length'] = query['end'] - query['start']
        return query

    def flanking_query(self, flank):
        query = self.flanking_reference_query(flank)
        if query is None:
            query = self.flanking_accession_query(flank)
        return query

    def to_fasta(self, flank, limit):
        sequence = self.sequence
        if flank is not None and flank > 0:
            flanking = self.flanking_query(flank)
            if flanking is not None:
                sequence = flanking['sequence']
        return to_fasta(self.id, sequence.nucleotide, None, limit)

    def validate(self):
        self.validate_in_accesion()
        self.validate_in_reference()

    def validate_in_accesion(self):
        accession = self.pipeline.resolver.accession_fetch(self.accession)
        if accession is not None:
            flanking = self.flanking_accession_query(0)
            if flanking is not None:
                if not self.sequence.nucleotide:
                    # if no sequence defined, assign sequence from accession 
                    self.sequence.nucleotide = flanking['sequence'].nucleotide
                    self.log.info('sequence assigned to artifact %s from accession %s:%d:%d', self.id, accession.id, self.accession_start, self.accession_end)
                else:
                    if flanking['sequence'].nucleotide != self.sequence.nucleotide:
                        self.log.info('artifact sequence in %s does not match accession %s:%d:%d', self.id, accession.id, self.accession_start, self.accession_end)
                        # self.search_in_accession()
                    else:
                        self.log.debug('artifact sequence for %s matched to accession %s:%d:%d', self.id, accession.id, self.accession_start, self.accession_end)

                    if self.strain is None and accession.strain is not None:
                        self.head['strain'] = accession.strain
                        self.log.info('strain %s assigned to artifact %s from %s', self.strain, self.id, self.accession)
            else:
                self.log.error('coordinates %d:%d invalid for accession %s', self.accession_start, self.accession_end, self.accession)
        else:
            self.log.error('accession %s could not be located', self.accession)

    def validate_in_reference(self):
        if self.aligned:
            flanking = self.flanking_reference_query(0)
            if flanking is not None:
                if not self.sequence.nucleotide:
                    # if no sequence defined, assign sequence from the reference 
                    self.sequence.nucleotide = flanking['sequence'].nucleotide
                    self.log.info('sequence assigned to artifact %s from reference %d:%d', self.id, self.reference_start, self.reference_end)
                else:
                    if flanking['sequence'].nucleotide != self.sequence.nucleotide:
                        self.log.info('artifact sequence in %s does not match reference %d:%d', self.id, self.reference_start, self.reference_end)
                    else:
                        self.log.debug('artifact sequence for %s matched to reference %d:%d', self.id, self.reference_start, self.reference_end)
            else:
                self.log.error('coordinates %d:%d invalid for reference', self.reference_start, self.reference_end)

    def search_in_accession(self):
        accession = self.pipeline.resolver.accession_fetch(self.accession)
        positions = []
        start = accession.sequence.nucleotide.find(self.sequence.nucleotide)
        while start >= 0 and start < accession.sequence.length:
            end = start + self.sequence.length
            positions.append([start, end])
            start = accession.sequence.nucleotide.find(self.sequence.nucleotide, start + 1)

        if len(positions) < 1:
            self.log.error('%s not found in %s', self.id, self.accession)
        elif len(positions) == 1:
            position = positions[0]
            self.log.info('new position found for %s in %s:%d:%d delta is %d:%d',
                self.id,
                self.accession,
                position[0],
                position[1],
                self.accession_start - position[0],
                self.accession_end - position[1])

            self.body['accession start'] = position[0]
            self.body['accession end'] = position[1]
        else:
            self.log.error('%s found %d times in %s', self.id, len(positions), self.accession)
            self.log.debug(str(positions))

    def align_to_reference(self, start, end, strand):
        self.log.info('gene %s aligned to %d:%d on the %s strand', self.id, start, end, 'plus' if strand else 'minus')
        self.body['reference start'] = start
        self.body['reference end'] = end
        self.body['reference strand'] = strand
        self.head['aligned'] = True

    def assign_reference_alignment(self, hit):
        self.align_to_reference(hit['target start'], hit['target end'], hit['query strand'])
        self.body['alignment'] = hit
        self.body['reference file sha1'] = hit['reference file sha1']
        self.body['reference record sha1'] = hit['reference record sha1']

class Gene(Artifact):
    def __init__(self, pipeline, node=None):
        Artifact.__init__(self, pipeline, node)
        self.log = logging.getLogger('Gene')
        for k in [ 'confirmed',  'identified' ]:
            if k not in self.head:
                self.head[k] = True

        for name in [
            '3 heptamer',
            '3 nonamer',
            '3 spacer',
            '5 heptamer',
            '5 nonamer',
            '5 spacer',
            '3 gap',
            '5 gap',
            'preamble',
            'fr1',
            'cdr1',
            'fr2',
            'cdr2',
            'fr3',
            'cdr3'
        ]:
            if name in self.body:
                artifact = self.body[name]
                if 'sequence' in artifact and isinstance(artifact['sequence'], dict):
                    artifact['sequence'] = Sequence(self.pipeline, artifact['sequence'])

                if 'start' in artifact and 'end' not in artifact and 'sequence' in artifact:
                    artifact['end'] = artifact['start'] + artifact['sequence'].length

    @property
    def row(self):
        return '| ' + ' | '.join(['{:<13}', '{:<12}', '{:<1}', '{:<3}', '{:<12}']).format(
            str(self.body['gene']), 
            str(self.body['family']), 
            str(self.body['functionality']), 
            str(self.body['length']), 
            str(self.body['reference start']), 
        )

    @property
    def framed(self):
        return self.head['framed']

    @property
    def functionality(self):
        return self.head['functionality']

    def validate_artifact(self):
        self.validate_artifact_in_gene()
        self.validate_artifact_in_reference()

    def validate_artifact_in_reference(self):
        for name in [
            '3 heptamer',
            '3 nonamer',
            '3 spacer',
            '5 heptamer',
            '5 nonamer',
            '5 spacer',
            '3 gap',
            '5 gap',
            'preamble',
            'fr1',
            'cdr1',
            'fr2',
            'cdr2',
            'fr3',
            'cdr3'
        ]:
            if name in self.body:
                node = self.body[name]
                if ('reference strand' in node and \
                    'reference start' in node and \
                    'reference end' in node):
                    artifact = Artifact(self.pipeline, { 'head':{ 'aligned':True, 'id': '{} {}'.format(name, self.id) }, 'body':node })
                    artifact.validate_in_reference()

    def validate_artifact_in_gene(self):
        for artifact in [
            'preamble',
            'fr1',
            'cdr1',
            'fr2',
            'cdr2',
            'fr3',
            'cdr3'
        ]:
            if artifact in self.body:
                start = self.body[artifact]['start']
                end = None if 'end' not in self.body[artifact] else self.body[artifact]['end']
                fragment = self.sequence.crop(start, end)
                if fragment:
                    self.body[artifact]['sequence'] = fragment
                    if artifact not in [ 'preamble', 'cdr3' ] and fragment.read_frame != 0:
                        self.log.debug('artifact %s for %s is out of frame %s', artifact, self.id, fragment.read_frame)
                else:
                    del self.body[artifact]
                    self.log.info('dropping empty %s for %s', artifact, self.id)

    def check_rss(self, flank, distance):
        for artifact in [
            '3 heptamer',
            '3 nonamer',
            '3 spacer',
            '5 heptamer',
            '5 nonamer',
            '5 spacer',
            '3 gap',
            '5 gap',
        ]:
            if artifact in self.body and self.body[artifact]['source'] in ['implied', 'pattern']:
                del self.body[artifact]

        flanking = self.flanking_reference_query(flank)
        if flanking:
            if self.region == 'vh':
                self._check_rss_3(flanking, 'vh', distance)

            elif self.region == 'dh':
                self._check_rss_5(flanking, 'dh', distance)
                self._check_rss_3(flanking, 'dh', distance)

            elif self.region == 'jh':
                self._check_rss_5(flanking, 'jh', distance)

    def html(self):
        buffer = []
        buffer.append('<div class="section">')
        buffer.append('<div class="genename">{} {}</div>'.format(self.id, self.functionality))

        self._codon_html(buffer)        
        if self.region == 'vh':
            self._nucleotide_html(buffer)
            self._html_3(buffer, 'vh')

        elif self.region == 'dh':
            self._html_5(buffer, 'dh')
            self._nucleotide_html(buffer)
            self._html_3(buffer, 'dh')

        elif self.region == 'jh':
            self._html_5(buffer, 'jh')
            self._nucleotide_html(buffer)
        buffer.append('</div>')
        print(''.join(buffer))

    def _codon_html(self, buffer):
        framed_artifact = True
        for artifact in [
            'fr1',
            'cdr1',
            'fr2',
            'cdr2',
            'fr3',
            'cdr3'
        ]:
            if artifact in self.body and self.body[artifact]['sequence'].read_frame != 0:
                framed_artifact = False

        if self.framed:
            buffer.append('<div class="frame codon">')
            if self.region == 'vh' and framed_artifact:
                for artifact in [
                    'preamble',
                    'fr1',
                    'cdr1',
                    'fr2',
                    'cdr2',
                    'fr3',
                    'cdr3'
                ]:
                    if artifact in self.body:
                        buffer.append('<span class="{}">{}</span>'.format(artifact, self.body[artifact]['sequence'].codon))
            else:
                buffer.append(self.sequence.codon)
            buffer.append('</div>')
        else:
            for frame in [0,1,2]:
                if self.sequence.read_frame == frame:
                    buffer.append('<div class="frame codon">')
                    buffer.append(self.sequence.codon_at_frame(frame))
                    buffer.append('</div>')
                else:
                    buffer.append('<div class="codon">')
                    buffer.append(self.sequence.codon_at_frame(frame))
                    buffer.append('</div>')

    def _nucleotide_html(self, buffer):
        buffer.append('<div class="gene">')
        if self.region == 'vh':
            for artifact in [
                'preamble',
                'fr1',
                'cdr1',
                'fr2',
                'cdr2',
                'fr3',
                'cdr3'
            ]:
                if artifact in self.body:
                    buffer.append('<span class="{}">{}</span>'.format(artifact, self.body[artifact]['sequence'].nucleotide))
        else:
            buffer.append(self.sequence.nucleotide)
        buffer.append('</div>')

    def _check_rss_3(self, flanking, region, distance):
        position = self._lookup_artifact(flanking, region, True, 0 , 'heptamer')
        position = self._lookup_artifact(flanking, region, True, position , 'nonamer')
        self._lookup_spacer(flanking, True)
        self._align_to_rss(flanking, True, distance)
        self._check_gap(flanking, True)

    def _check_rss_5(self, flanking, region, distance):
        position = self._lookup_artifact(flanking, region, False, 0, 'heptamer')
        position = self._lookup_artifact(flanking, region, False, position , 'nonamer')
        self._lookup_spacer(flanking, False)
        self._align_to_rss(flanking, False, distance)
        self._check_gap(flanking, False)

    def _check_gap(self, flanking, orientation):
        gap = '{} gap'.format('3' if orientation else '5')
        heptamer = '{} heptamer'.format('3' if orientation else '5')
        nonamer = '{} nonamer'.format('3' if orientation else '5')

        if gap in self.body: del self.body[gap]
        offset = 0
        if heptamer in self.body:
            if orientation == self.body['reference strand']:
                offset = self.body[heptamer]['reference start'] - self.body['reference end']
            else:
                offset = self.body['reference start'] - self.body[heptamer]['reference end'] 

        elif nonamer in self.body:
            if orientation == self.body['reference strand']:
                offset = self.body[nonamer]['reference start'] - self.body['reference end']
            else:
                offset = self.body['reference start'] - self.body[nonamer]['reference end'] 

        if offset > 0:
            artifact = {
                'reference strand': flanking['reference strand'],
                'source': 'implied',
            }

            if orientation == self.body['reference strand']:
                artifact['reference start'] = self.body['reference end']
                artifact['reference end'] = self.body['reference end'] + offset
            else:
                artifact['reference start'] = self.body['reference start'] - offset
                artifact['reference end'] = self.body['reference start']

            if flanking['reference strand']:
                start = artifact['reference start'] - flanking['flank start']
                end = artifact['reference end'] - flanking['flank start']
            else:
                start = flanking['flank end'] - artifact['reference end']
                end = flanking['flank end'] - artifact['reference start']

            artifact['sequence'] = flanking['sequence'].crop(start, end)
            self.body[gap] = artifact
            self.log.info('annotating %s with %dbp %s %s', self.id, artifact['sequence'].length, gap, artifact['sequence'].nucleotide)

    def _align_to_rss(self, flanking, orientation, distance):
        name = '{} heptamer'.format('3' if orientation else '5')
        offset = 0
        if name in self.body:
            if orientation == self.body['reference strand']:
                offset = self.body[name]['reference start'] - self.body['reference end']
            else:
                offset = self.body['reference start'] - self.body[name]['reference end'] 

            if offset > 0:
                if distance is None or offset <= distance:
                    if orientation == self.body['reference strand']:
                        self.body['reference end'] += offset
                    else:
                        self.body['reference start'] -= offset

                    if orientation == self.accession_strand:
                        self.body['accession end'] += offset
                    else:
                        self.body['accession start'] -= offset

                    self.log.info('aligning %s gene with %s by %dbp', self.id, name, offset)
                    self.body['length'] = self.accession_end - self.accession_start
                    self.sequence.nucleotide = None
                else:
                    self.log.info('not aligning %s gene with %s because distance is %dbp which is further than %dbp', self.id, name, offset, distance)

    def _lookup_artifact(self, flanking, region, orientation, offset, type):
        position = 0
        name = '{} {}'.format('3' if orientation else '5', type)
        pattern = self.configuration.rss[region][name]
        if name not in self.body:
            if orientation:
                # Look in the 3 prime flanking sequence
                interval = flanking['sequence'].crop(flanking['end'] + offset, flanking['flank length'])

            else:
                # Look for the reverse complement in the reverse complement of the 5 prime flanking sequence
                interval = flanking['sequence'].crop(0, flanking['start'] - offset).reversed

            match = pattern.search(interval.nucleotide)
            if match:
                artifact = {
                    'reference strand': flanking['reference strand'],
                    'source': 'pattern',
                    'sequence': Sequence(self.pipeline, {
                        'nucleotide': match.group(0),
                        'read frame': 0,
                        'strand': True
                    }),
                }
                if not orientation: artifact['sequence'] = artifact['sequence'].reversed

                if orientation == flanking['reference strand']:
                    artifact['reference start'] = flanking['reference end'] + match.start() + offset
                    artifact['reference end'] = flanking['reference end'] + match.end() + offset
                    gap = artifact['reference start'] - flanking['reference end']
                else:
                    artifact['reference start'] = flanking['reference start'] - match.end() - offset
                    artifact['reference end'] = flanking['reference start'] - match.start() - offset
                    gap = flanking['reference start'] - artifact['reference end']

                artifact['length'] = artifact['reference end'] - artifact['reference start']
                self.body[name] = artifact
                self.log.info('%s %s found %dbp from the %s gene', name, artifact['sequence'].nucleotide, gap, self.id)
                position = artifact['length'] + gap
        else:
            artifact = self.body[name]
            if orientation == flanking['reference strand']:
                artifact['reference start'] = flanking['reference end'] + artifact['relative start']
                artifact['reference end'] = flanking['reference end'] + artifact['relative end']
                gap = artifact['reference start'] - flanking['reference end']
            else:
                artifact['reference start'] = flanking['reference start'] - artifact['relative end']
                artifact['reference end'] = flanking['reference start'] - artifact['relative start']
                gap = flanking['reference start'] - artifact['reference end']
            artifact['length'] = artifact['reference end'] - artifact['reference start']
            self.log.info('%s %s aligned %dbp from the %s gene', name, artifact['sequence'].nucleotide, gap, self.id)
            position = artifact['length'] + gap
        return position

    def _lookup_spacer(self, flanking, orientation):
        name = '{} spacer'.format('3' if orientation else '5')
        heptamer_name = '{} heptamer'.format('3' if orientation else '5')
        nonamer_name = '{} nonamer'.format('3' if orientation else '5')
        if  name not in self.body and heptamer_name in self.body and nonamer_name in self.body:
            heptamer = self.body[heptamer_name]
            nonamer = self.body[nonamer_name]
            artifact = {
                'source': 'implied',
                'reference strand': flanking['reference strand'],
            }
            if flanking['reference strand'] == orientation:
                artifact['reference start'] = heptamer['reference end']
                artifact['reference end'] = nonamer['reference start']
            else:
                artifact['reference start'] = nonamer['reference end']
                artifact['reference end'] = heptamer['reference start']

            if flanking['reference strand']:
                start = artifact['reference start'] - flanking['flank start']
                end = artifact['reference end'] - flanking['flank start']
            else:
                start = flanking['flank end'] - artifact['reference end']
                end = flanking['flank end'] - artifact['reference start']

            artifact['sequence'] = Sequence(self.pipeline, {
                'nucleotide': flanking['sequence'].nucleotide[start:end],
                'strand': True,
                'read frame': 0
            })
            artifact['length'] = artifact['reference end'] - artifact['reference start']
            self.body[name] = artifact
            self.log.info('%s %s of length %dbp found for %s gene', name, artifact['sequence'].nucleotide, artifact['sequence'].length, self.id)

    def _html_3(self, buffer, region):
        if ('3 heptamer' in self.body or '3 nonamer' in self.body) and '3 gap' in self.body:
            buffer.append('<div class="gap">{}</div>'.format(self.body['3 gap']['sequence'].nucleotide))

        if '3 heptamer' in self.body:
            buffer.append('<div class="heptamer {}">{}</div>'.format(self.body['3 heptamer']['source'], self.body['3 heptamer']['sequence'].nucleotide)) 

        if '3 spacer' in self.body:
            buffer.append('<div class="spacer {}">{}</div>'.format(self.body['3 spacer']['source'], self.body['3 spacer']['sequence'].nucleotide))

        if '3 nonamer' in self.body:
            buffer.append('<div class="nonamer {}">{}</div>'.format(self.body['3 nonamer']['source'], self.body['3 nonamer']['sequence'].nucleotide))

    def _html_5(self, buffer, region):
        if '5 nonamer' in self.body:
            buffer.append('<div class="nonamer {}">{}</div>'.format(self.body['5 nonamer']['source'], self.body['5 nonamer']['sequence'].nucleotide))

        if '5 spacer' in self.body:
            buffer.append('<div class="spacer {}">{}</div>'.format(self.body['5 spacer']['source'], self.body['5 spacer']['sequence'].nucleotide))

        if '5 heptamer' in self.body:
            buffer.append('<div class="heptamer {}">{}</div>'.format(self.body['5 heptamer']['source'], self.body['5 heptamer']['sequence'].nucleotide)) 

        if ('5 heptamer' in self.body or '5 nonamer' in self.body) and '5 gap' in self.body:
            buffer.append('<div class="gap">{}</div>'.format(self.body['5 gap']['sequence'].nucleotide))

    def assign_reference_alignment(self, hit):
        Artifact.assign_reference_alignment(self, hit)
        for name in [
            'preamble',
            'fr1',
            'cdr1',
            'fr2',
            'cdr2',
            'fr3',
            'cdr3'
        ]:
            if name in self.body:
                artifact = self.body[name]
                artifact['reference file sha1'] = self.body['reference file sha1']
                artifact['reference record sha1'] = self.body['reference record sha1']

                if artifact['sequence'].strand == self.sequence.strand:
                    artifact['reference strand'] = self.reference_strand
                else:
                    artifact['reference strand'] = not self.reference_strand

                if artifact['reference strand']:
                    artifact['reference start'] = self.reference_start + artifact['start']
                    artifact['reference end'] = self.reference_start + artifact['end']
                else:
                    artifact['reference start'] = self.reference_start + self.length - artifact['end']
                    artifact['reference end'] = self.reference_start + self.length - artifact['start']

class Hit(object):
    def __init__(self, sample, node=None):
        self.log = logging.getLogger('Hit')
        self.sample = sample
        self._operation = None
        self.node = {
            'valid': True,
            'picked': False,

            'uuid': None,
            'type': None,
            'region': None,
            'FLAG': 0,
            'CIGAR': None,
            'MAPQ': None,
            'query start': None,
            'query end': None,
            'subject id': None,
            'subject start': None,
            'subject end': None,
            'query': None,
            'codon mismatch': None,
            'charge': None,
            'weight': None,
            'in frame': None,

            'strain': None,
            'gene': None,
            'family': None,
            'allele': None,
            'framed': None,
            'functionality': None,
            'subject': None,
            'subject strand': None,
            'reference start': None,
            'reference end': None,
        }
        self.node = merge(self.node, node)
        if self.node['uuid'] is None:
            self.node['uuid'] = str(uuid.uuid4())

    @property
    def pipeline(self):
        return self.sample.pipeline

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def document(self):
        document = to_document(self.node)
        return document

    @property
    def uuid(self):
        return self.node['uuid']

    @property
    def strain(self):
        if self.node['strain'] is None:
            self._load_subject()
        return self.node['strain']

    @property
    def gene(self):
        if self.node['gene'] is None:
            self._load_subject()
        return self.node['gene']

    @property
    def family(self):
        if self.node['family'] is None:
            self._load_subject()
        return self.node['family']

    @property
    def allele(self):
        if self.node['allele'] is None:
            self._load_subject()
        return self.node['allele']

    @property
    def framed(self):
        if self.node['framed'] is None:
            self._load_subject()
        return self.node['framed']

    @property
    def functionality(self):
        if self.node['functionality'] is None:
            self._load_subject()
        return self.node['functionality']

    @property
    def subject(self):
        if self.node['subject'] is None:
            self._load_subject()
        return self.node['subject']

    @property
    def query_strand(self):
        return self.subject_strand != self.reverse_complemented

    @property
    def subject_strand(self):
        if self.node['subject strand'] is None:
            self._load_subject()
        return self.node['subject strand']

    @property
    def reference_start(self):
        if self.node['reference start'] is None:
            self._load_subject()
        return self.node['reference start']

    @property
    def reference_end(self):
        if self.node['reference end'] is None:
            self._load_subject()
        return self.node['reference end']

    @property
    def FLAG(self):
        return self.node['FLAG']

    @property
    def reverse_complemented(self):
        return self.FLAG & 0x10 == 0x10

    @property
    def unmapped(self):
        return self.FLAG & 0x4 == 0x4

    @property
    def region(self):
        return self.node['region']

    @region.setter
    def region(self, value):
        self.node['region'] = value

    @property
    def subject_id(self):
        return self.node['subject id']

    @property
    def valid(self):
        return self.node['valid']

    @property
    def picked(self):
        return self.node['picked']

    @picked.setter
    def picked(self, value):
        self.node['picked'] = value

    @property
    def sufficient(self):
        return self.valid

    @property
    def type(self):
        return self.node['type']

    @type.setter
    def type(self, value):
        self.node['type'] = value

    @property
    def palindrome(self):
        return self.node['palindrome']

    @property
    def primary(self):
        return self.node['primary']

    @primary.setter
    def primary(self, value):
        self.node['primary'] = value

    @property
    def in_frame(self):
        if self.node['in frame'] is None:
            self.node['in frame'] = (self.query.read_frame == self.subject.read_frame)
        return self.node['in frame']

    @property
    def gapped(self):
        if 'gapped' not in self.node:
            self.node['gapped'] = ( 'D' in self.CIGAR or 'I' in self.CIGAR )
        return self.node['gapped']

    @property
    def score(self):
        return self.node['score']

    @score.setter
    def score(self, value):
        self.node['score'] = value

    @property
    def subject_start(self):
        return self.node['subject start']

    @property
    def subject_end(self):
        return self.node['subject end']

    @property
    def query_start(self):
        return self.node['query start']

    @property
    def query_end(self):
        return self.node['query end']

    @property
    def query(self):
        if self.node['query'] is None:
            # when the 0x10 BAM flag is on the alignments refer to the reverse complement of the read
            if self.reverse_complemented:
                self.node['query'] = self.sample.sequence.reversed.crop(self.query_start, self.query_end)
            else:
                self.node['query'] = self.sample.sequence.crop(self.query_start, self.query_end)
        return self.node['query']

    @property
    def chew_3_prime(self):
        if '3 chew' not in self.node:
            if  self.valid and \
                self.type == 'gene segment' and \
                (self.region == 'vh' or self.region == 'dh'):
                gene = self.pipeline.resolver.gene_fetch(self.subject_id)
                if gene:
                    if self.subject_end < gene.sequence.length:
                        self.node['3 chew'] = gene.sequence.crop(self.subject_end, gene.sequence.length)
        return None if '3 chew' not in self.node else self.node['3 chew']

    @property
    def chew_5_prime(self):
        if '5 chew' not in self.node:
            if  self.valid and \
                self.type == 'gene segment' and \
                (self.region == 'jh' or self.region == 'dh'):
                gene = self.pipeline.resolver.gene_fetch(self.subject_id)
                if gene:
                    if self.subject_start > 0:
                        self.node['5 chew'] = gene.sequence.crop(0, self.subject_start)
        return None if '5 chew' not in self.node else self.node['5 chew']

    @property
    def length(self):
        return self.query_end - self.query_start

    @property
    def MAPQ(self):
        return self.node['MAPQ']

    @property
    def CIGAR(self):
        return self.node['CIGAR']

    @property
    def operation(self):
        if self._operation is None and self.CIGAR:
            if self.configuration.expression['cigar'].fullmatch(self.CIGAR):
                self._operation = []
                if self.CIGAR != '*':
                    end = len(self.CIGAR) - 1
                    p = 0
                    o = 0
                    while p < end:
                        if self.CIGAR[p] in 'MIDNSHPX=':
                            self._operation.append((self.CIGAR[p], int(self.CIGAR[o:p])))
                            o = p + 1
                        p += 1
                    self._operation.append((self.CIGAR[-1:], int(self.CIGAR[o:-1])))
            else:
                raise ValueError('expression {} is not a valid CIGAR string'.format(self.CIGAR))
        return self._operation

    @property
    def TLEN(self):
        return self.node['TLEN']

    @property
    def codon_mismatch(self):
        if self.node['codon mismatch'] is None and self.subject.codon and self.query.codon:
            self.node['codon mismatch'] = 0
            for s,q in zip(self.subject.codon, self.query.codon):
                if s != q:
                    self.node['codon mismatch'] += 1
        return self.node['codon mismatch']

    @property
    def charge(self):
        if self.node['charge'] is None:
            if self.query.codon:
                self.node['charge'] = 0
                for acid in self.query.codon:
                    amino = self.configuration.enumeration['iupac amino acid notation']
                    if 'charge' in amino[acid]:
                        self.node['charge'] += amino[acid]['charge']
        return self.node['charge']

    @property
    def weight(self):
        if self.node['weight'] is None:
            if self.query.codon:
                self.node['weight'] = 0
                for acid in self.query.codon:
                    amino = self.configuration.enumeration['iupac amino acid notation']
                    if 'weight' in amino[acid]:
                        self.node['weight'] += amino[acid]['weight']
        return self.node['weight']

    @classmethod
    def from_bwa(cls, block, value, region, type='gene segment'):
        instance = None
        if block is not None:
            if value is not None:
                # 0  QNAME   string     Query template NAME
                # 1  FLAG    int        bitwise FLAG
                #    0x1     template having multiple segments in sequencing
                #    0x2     each segment properly aligned according to the aligner
                #    0x4     segment unmapped
                #    0x8     next segment in the template unmapped
                #    0x10    SEQ being reverse complemented
                #    0x20    SEQ of the next segment in the template being reverse complemented
                #    0x40    the first segment in the template
                #    0x80    the last segment in the template
                #    0x100   secondary alignment
                #    0x200   not passing filters, such as platform/vendor quality controls
                #    0x400   PCR or optical duplicate
                #    0x800   supplementary alignment
                #
                # 2  RNAME   string     Reference sequence NAME
                # 3  POS     int        1-based leftmost mapping POSition
                # 4  MAPQ    int        MAPping Quality
                # 5  CIGAR   string     CIGAR string
                # 6  RNEXT   string     Reference name of the mate/next read
                # 7  PNEXT   int        Position of the mate/next read
                # 8  TLEN    int        observed Template LENgth
                # 9  SEQ     string     segment SEQuence
                # 10 QUAL    string     Phred QUALity+33
                #
                # MD    Z   String for mismatching positions. [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)
                # AS    i   Alignment score generated by aligner
                # NM    i   Edit distance to the reference, including ambiguous bases but excluding clipping
                # XS        Suboptimal alignment score
                match = block.configuration.expression['sam record'].search(value)
                if match:
                    record = match.groupdict()
                    record['region'] = region
                    record['type'] = type
                    for field in [
                        'FLAG',
                        'POS',
                        'MAPQ',
                        'PNEXT',
                        'TLEN'
                    ]:
                        record[field] = int(record[field])

                    if record['QNAME'] in block.lookup:
                        if record['FLAG'] & 0x4 == 0:
                            sample = block.lookup[record['QNAME']]
                            record['subject id'] = record['RNAME']
                            record['subject start'] = record['POS'] - 1
                            if 'OPT' in record:
                                optional = record['OPT'].strip('\t').split('\t')
                                for o in optional:
                                    m = block.configuration.expression['sam optional field'].search(o)
                                    if m:
                                        o = m.groupdict()
                                        if block.configuration.expression['sam optional field value'][o['TYPE']].match(o['VALUE']):
                                            if o['TYPE'] is 'i':
                                                o['VALUE'] = int(o['VALUE'])
                                            elif o['TYPE'] is 'f':
                                                o['VALUE'] = float(o['VALUE'])
                                            record[o['TAG']] = o['VALUE']
                                        else:
                                            block.log.error('ignoring invalid %s optional field %s', o['TYPE'], o['VALUE'])

                            if 'AS' in record:
                                record['score'] = record['AS']

                            for field in [
                                'QUAL',
                                'SEQ',
                                'QNAME',
                                'RNAME',
                                'RNEXT',
                                'PNEXT',
                                'OPT',
                                'POS'
                            ]:
                                if field in record:
                                    del record[field]

                            record['subject end'] = record['subject start']
                            record['query start'] = 0
                            record['query end'] = 0
                            instance = Hit(sample, record)
                            if instance.operation:
                                inside = False
                                for o in instance.operation:
                                    if not inside:
                                        if o[0] in ['S', 'H']:
                                            instance.node['query start'] += o[1]
                                        else:
                                            inside = True
                                            instance.node['query end'] = instance.node['query start']
                                    if inside:
                                        if o[0] in ['M', '=', 'X']:
                                            instance.node['query end'] += o[1]
                                            instance.node['subject end'] += o[1]
                                        elif o[0] is 'I':
                                            pass
                                        elif o[0] is 'D':
                                            instance.node['subject end'] += o[1]
                                        elif o[0] in ['S', 'H']:
                                            break

                            if not instance.query_start < instance.query_end:
                                block.log.error('discarding zero length hit %s', value)
                                instance = None
                        else:
                            block.log.debug('discarding unmapped hit for sample %s', record['QNAME'])
                    else:
                        block.log.error('could not locate sample %s in block', record['QNAME'])
                else:
                    block.log.error('invalid sam syntax %s', value)
            else:
                block.log.error('sam record cannot be null')
        else:
            block.log.error('block cannot be null')
        return instance

    @classmethod
    def from_igblast(cls, block, value, qname):
        instance = None
        if block is not None:
            if value is not None:
                match = self.configuration.expression['igblast hit'].search(value)
                if match:
                    record = dict(((k.replace('_', ' '),v) for k,v in match.groupdict().items() if k and v))
                    record['QNAME'] = qname
                    if record['subject strand'] == 'plus':
                        record['subject strand'] = True
                    elif record['subject strand'] == 'minus':
                        record['subject strand'] = False
                    else:
                        self.log.warning('could not parse value %s for as a strand setting to plus', record['subject strand'])
                        record['subject strand'] = True
                    for k in [
                        'identical',
                        'bit score',
                        'evalue',
                    ]:
                        try:
                            record[k] = float(record[k])
                        except ValueError:
                            self.log.warning('could not parse value %s for %s as int', record[k], k)
                            record[k] = None

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
                        try: 
                            record[k] = int(record[k])
                        except ValueError as e:
                            self.log.warning('could not parse value %s for %s as int', record[k], k)
                            record[k] = None

                    # fix the region for heavy chain
                    if 'region' in record:
                        record['region'] = record['region'].lower() + 'h'

                    # switch from % to fraction
                    if 'identical' in record:
                        record['identical'] /= 100.0

                    if record['QNAME'] in block.lookup:
                        sample = block.lookup[record['QNAME']]
                        instance = Hit(sample, record)
            else:
                block.log.error('igblast record cannot be null')
        else:
            block.log.error('block cannot be null')
        return instance

    @classmethod
    def clone(cls, hit):
        instance = None
        if hit is not None:
            node = {}
            for key, value in hit.node.items():
                if key not in [
                    'uuid',
                    'query',
                    'subject',
                    'reference start',
                    'reference end',
                    'subject strand',
                ]:
                    node[key] = deepcopy(value)
                else:
                    node[key] = None

            instance = Hit(hit.sample, node)
        return instance

    def invalidate(self, message):
        self.node['valid'] = False
        if 'comment' not in self.node:
            self.node['comment'] = []
        self.node['comment'].append(message)
        self.log.info('%dbp %s hit to %s on %s : %s', self.length, self.region, self.subject_id, self.sample.id, message)

    def trim(self, start, end):
        start = max(start, self.query_start)
        end = min(end, self.query_end)
        if start < end:
            if start > self.query_start:
                overlap = start - self.query_start
                if overlap > 0:
                    self.node['query start'] = start
                    self.node['subject start'] += overlap
                    self.score -= overlap * self.configuration.selected['region']['dh']['overlap penalty factor']
                    self.reset()

            if end < self.query_end:
                overlap = self.query_end - end
                if overlap > 0:
                    self.node['query end'] = end
                    self.node['subject end'] -= overlap
                    self.score -= overlap * self.configuration.selected['region']['dh']['overlap penalty factor']
                    self.reset()
        else:
            self.invalidate('hit trimmed to length zero')
        return self.valid

    def nonoverlapping(self, start, end):
        result = 0
        start = max(start, self.query_start)
        end = min(end, self.query_end)
        if start < end:
            overlap = 0
            if start > self.query_start:
                overlap += start - self.query_start
            if end < self.query_end:
                overlap += self.query_end - end
            result = self.query_end - self.query_start - overlap
        return result

    def reset(self):
        self.node['valid'] = True
        self.node['picked'] = False
        for term in [
            'charge',
            'weight',
            'strain',
            'gene',
            'family',
            'allele',
            'functionality',
            'framed',
            'in frame',
            'query',
            'subject',
            'subject strand'
            'reference start',
            'reference end',
        ]:
            if term in self.node:
                self.node[term] = None

    def _load_subject(self):
        if self.type == 'gene segment':
            gene = self.pipeline.resolver.gene_fetch(self.subject_id)
            if gene:
                self.node['strain'] = gene.strain
                self.node['gene'] = gene.body['gene']
                self.node['family'] = gene.body['family']
                self.node['allele'] = gene.body['allele']
                self.node['framed'] = gene.framed
                self.node['functionality'] = gene.functionality
                self.node['subject'] = gene.sequence.crop(self.subject_start, self.subject_end)

                # this my need to be different for gapped alignments
                self.query.read_frame = self.node['subject'].read_frame

                if gene.aligned:
                    if gene.reference_strand:
                        self.node['reference start'] = gene.reference_start + self.subject_start
                        self.node['reference end'] = gene.reference_start + self.subject_end
                    else:
                        self.node['reference start'] = gene.reference_start + gene.length - self.subject_end
                        self.node['reference end'] = gene.reference_start + gene.length - self.subject_start
                    self.node['subject strand'] = gene.reference_strand

class Sample(object):
    def __init__(self, pipeline, node=None, id=None, library=None):
        self.log = logging.getLogger('Sample')
        self.pipeline = pipeline
        self.node = node
        self.index = {}
        self._primary = None
        self._observation = None
        
        if self.node is None:
            self.node = {
                'head': {
                    'id': None,
                    'FLAG': 0,
                    'library': library,
                    'abundance': 1,
                    'valid': True,
                }
            }

            if id:
                fragmented = id.split()
                self.head['id'] = fragmented[0]
                if len(fragmented) > 1:
                    self.head['id comment'] = ' '.join(fragmented[1:])
                    try:
                        self.head['abundance'] = int(fragmented[1].split(':')[0])
                    except ValueError:
                        pass
            self._reset()

        if not self.id:
            raise InvalidSampleError('sample must have an id')
            
        if not self.key:
            raise InvalidSampleError('sample must have a valid key')
            
        if not self.library:
            raise InvalidSampleError('sample must have a library')
            
        for sequence in ['sequence', 'effective']:
            if sequence in self.body:
                if isinstance(self.body[sequence], dict):
                    self.body[sequence] = Sequence(self.pipeline, self.body[sequence])
            else:
                self.body[sequence] = Sequence(self.pipeline)
            
        for hit in self.hit:
            self._load_hit(hit)

    def __str__(self):
        return self.id

    def add_bwa_hit(self, hit):
        if hit is not None and not hit.unmapped:
            gene = self.pipeline.resolver.gene_fetch(hit.subject_id)
            if gene is None:
                hit.invalidate('unknown gene')
            self._add_hit(hit)

    def add_igblast_hit(self, hit):
        if hit is not None:
            if 'subject id' in hit:
                gene = self.pipeline.resolver.gene_fetch(hit.subject_id)
                # switch the hit to a zero based coordinate system
                if gene:
                    hit.query_start -= 1
                    if hit.node['subject strand'] == gene.sequence.strand:
                        hit.subject_start -= 1
                    else:
                        hit.subject_start = gene.length - hit.subject_start
                        hit.subject_end = gene.length - hit.subject_end + 1
                else:
                    hit.invalidate('unknown gene')
            self._add_hit(hit)

    @property
    def residual(self):
        start = 0
        end = self.sequence.length
        if 'jh' in self.primary:
            jh = self.primary['jh']
            if jh.query_strand:
                if jh.query.strand != self.sequence.strand:
                    start = max(start, self.sequence.length - jh.query_start)
                else:
                    end = min(end, jh.query_start)
            else:
                if jh.query.strand != self.sequence.strand:
                    start = max(start, self.sequence.length - jh.query_start)
                else:
                    end = min(end, jh.query_start)

        if 'vh' in self.primary:
            vh = self.primary['vh']
            if vh.query_strand:
                if vh.query.strand != self.sequence.strand:
                    end = min(end, self.sequence.length - vh.query_end)
                else:
                    start = max(start, vh.query_end)
            else:
                if vh.query.strand != self.sequence.strand:
                    end = min(end, self.sequence.length - vh.query_end)
                else:
                    start = max(start, vh.query_end)

        return self.sequence.crop(start, end)

    def determine_region(self, region):
        if region == 'jh':
            self._pick_jh_region()
        elif region == 'vh':
            self._pick_vh_region()
        elif region == 'dh':
            self._pick_dh_region()

    def analyze(self):
        self._reset()
        self._load_gene()
        self._pick_jh_region()
        self._pick_vh_region()
        self._check_v_j_framing()
        self._pick_dh_region()
        self._identify_v_j_junction()
        self._identify_v_d_junction()
        self._identify_d_j_junction()
        self._identify_cdr3()
        self._identify_chewback()
        self._check_for_stop_codon()
        self._check_for_productivity()
        self._analyze_quality()

    def _add_hit(self, hit):
        self.hit.append(hit)
        self._load_hit(hit)
        for term in [
            'effective',
            'effective query start',
            'effective query end'
        ]:
            if term in self.body:
                del self.body[term]

    def _load_hit(self, hit):
        if hit.uuid not in self.index:
            self.index[hit.uuid] = hit

        for k in [
            'query',
            'subject',
            '3 chew',
            '5 chew'
        ]:
            if k in hit.node and isinstance(hit.node[k], dict):
                hit.node[k] = Sequence(self.pipeline, hit.node[k])

    def _reset(self):
        self.head['valid'] = True
        self.head['gapped'] = False
        self.head['framed'] = False
        self.head['premature'] = False
        self.head['in frame'] = False
        self.head['productive'] = False
        self.head['palindromic'] = False
        self.head['p count'] = 0
        self.head['n count'] = 0
        self.head['conserved cdr3 c'] = False
        self.head['conserved cdr3 w'] = False
        self.head['cdr3 length'] = None
        self.head['dh length'] = None
        self.head['jh length'] = None
        self.head['vh length'] = None
        self.head['average phred'] = None
        self.body['primary'] = {}

        for term in [
            'effective',
            'comment',
            'framed by',
            'effective query start',
            'effective query end'
        ]:
            if term in self.body:
                del self.body[term]

        # remove all implicit hits
        if 'hit' in self.body:
            self.body['hit'] = [ hit for hit in self.body['hit'] if hit.type == 'gene segment']

        for hit in self.hit:
            hit.reset()

        self.head['gapped'] = any([ hit.gapped for hit in self.hit ])

    def _load_gene(self):
        for hit in self.hit:
            if hit.valid and hit.type == 'gene segment':
                if hit.subject is None:
                    hit.invalidate('unknown gene')

    def _pick_jh_region(self):
        return self._pick_region('jh', self.hit)

    def _pick_vh_region(self):
        return self._pick_region('vh', self.hit)

    def _pick_dh_region(self):
        jh = None if 'jh' not in self.primary else self.primary['jh']
        vh = None if 'vh' not in self.primary else self.primary['vh']

        if vh is not None and jh is not None:
            top = None
            search = { 'valid': [], 'map': {} }
            for hit in self.hit:
                if hit.region == 'dh':
                    search['map'][hit.uuid] = hit
                    if hit.valid:
                        search['valid'].append(hit)

            if len(search['valid']) > 0:
                search['trimmed'] = []
                for hit in search['valid']:
                    length = hit.nonoverlapping(vh.query_end, jh.query_start)
                    if length > 0:
                        if length < hit.length:
                            trimmed = Hit.clone(hit)
                            trimmed.trim(vh.query_end, jh.query_start)
                            trimmed.node['trimmed'] = True
                            self.add_bwa_hit(trimmed)
                            search['trimmed'].append(trimmed)
                        else:
                            search['trimmed'].append(hit)

                if len(search['trimmed']) > 0:
                    top = search['trimmed']

                self._pick_region('dh', top)

    def _pick_region(self, region, candidate):
        picked = False
        if candidate:
            top = None
            search = { 'valid': [] }
            for hit in candidate:
                if hit.region == region and hit.sufficient:
                    search['valid'].append(hit)

            if search['valid']:
                top = search['valid']
                top.sort(reverse=True, key=lambda x: x.score)
                top = [ hit for hit in top if hit.score == top[0].score ]

                if len(top) > 1:
                    search['framed'] = []
                    for hit in top:
                        if hit.framed and hit.codon_mismatch is not None:
                            search['framed'].append(hit)

                    if search['framed']:
                        top = search['framed']
                        top.sort(key=lambda x: x.codon_mismatch)
                        top = [ hit for hit in top if hit.codon_mismatch == top[0].codon_mismatch ]

                # if there are multiple matches and at least one is functional, keep only the functional
                if len(top) > 1:
                    search['functional'] = [ hit for hit in top if hit.functionality == 'F' ]
                    if search['functional']:
                        top = search['functional']
            if top:
                # prefer a functional gene for a primary F > O > P
                top.sort(key=lambda x: x.functionality)
                for hit in top:
                    self._pick_hit(hit)
                picked = True
        return picked

    def _pick_hit(self, hit):
        def pick_region_primary(hit):
            if hit.region not in self.body['primary']:
                hit.primary = True
                self._primary = None
                self.body['primary'][hit.region] = hit.uuid

        def pick_framing_hit(framing):
            if  not self.framed and \
                framing.framed and \
                framing.query_start is not None and \
                framing.subject is not None:

                # first hit to be picked that is framed sets the frame for the sample
                # that also means deciding if the entire sample is reverse complemented
                # any further hits will have to have the same orientation (except for dh region)
                self.head['framed'] = True
                self.body['framed by'] = framing.uuid
                self.reverse_complemented = framing.reverse_complemented
                self.sequence.read_frame = (framing.query_start + framing.subject.read_frame) % 3
                if self.reverse_complemented:
                    self.sequence.read_frame = (self.sequence.length - self.sequence.read_frame) % 3

                for hit in self.hit:
                    if hit.valid and hit.uuid != framing.uuid:
                        if self.reverse_complemented:
                            hit.query.read_frame = 2 - (hit.query_start - self.sequence.reversed.read_frame - 1) % 3
                        else:
                            hit.query.read_frame = 2 - (hit.query_start - self.sequence.read_frame - 1) % 3

                # self.log.debug('%s framed by %s', str(self), to_json_document(hit))

        if hit is not None and hit.query is not None:
            hit.picked = True
            pick_region_primary(hit)
            pick_framing_hit(hit)

    def _check_v_j_framing(self):
        jh = None if 'jh' not in self.primary else self.primary['jh']
        vh = None if 'vh' not in self.primary else self.primary['vh']
        if jh is not None and vh is not None:
            if jh.reverse_complemented == vh.reverse_complemented:
                if jh.in_frame and vh.in_frame:
                    self.head['in frame'] = True
            else:
                self.invalidate('vh and jh are on opposite strands')
        else:
            self.invalidate('could not establish a vh jh pair')

    def _identify_junction(self, left, right, name):
        def identify_palindrome(junction, left, right):
            if junction.query.length > 0:
                junction.node['palindrome'] = [ 'N' ] * junction.query.length
                i = 0
                while i < left.length and i < junction.query.length:
                    mirror = -(i + 1)
                    nucleotide = junction.query.nucleotide[i]
                    palindrome = left.nucleotide[mirror]
                    if self.configuration.enumeration['reverse complement'][nucleotide] == palindrome:
                        junction.node['palindrome'][i] = 'P'
                        i += 1
                    else:
                        break
                i = 0
                while i < right.length and i < junction.query.length:
                    mirror = -(i + 1)
                    nucleotide = right.nucleotide[i]
                    palindrome = junction.query.nucleotide[mirror]
                    if self.configuration.enumeration['reverse complement'][nucleotide] == palindrome:
                        junction.node['palindrome'][mirror] = 'P'
                        i += 1
                    else:
                        break
                junction.node['palindrome'] = ''.join(junction.node['palindrome'])
                p = junction.node['palindrome'].count('P')
                n = junction.node['palindrome'].count('N')
                self.head['p count'] += p
                self.head['n count'] += n
                if p > 0:
                    self.head['palindromic'] = True
                if n > 0:
                    if p > 0:
                        junction.node['palindrome ratio'] = float(n) / float(p)
                    else:
                        junction.node['palindrome ratio'] = 1.0
                else:
                    junction.node['palindrome ratio'] = 0.0

        if left is not None and right is not None and left.query_end < right.query_start:
            junction = Hit(self, {
                'FLAG': self.FLAG,
                'region': name,
                'subject id': name,
                'query start': left.query_end,
                'query end': right.query_start,
            })
            identify_palindrome(junction, left.query, right.query)
            self._add_hit(junction)
            self._pick_hit(junction)

    def _identify_v_d_junction(self):
        vh = None if 'vh' not in self.primary else self.primary['vh']
        dh = None if 'dh' not in self.primary else self.primary['dh']
        self._identify_junction(vh, dh, 'V-D')

    def _identify_d_j_junction(self):
        dh = None if 'dh' not in self.primary else self.primary['dh']
        jh = None if 'jh' not in self.primary else self.primary['jh']
        self._identify_junction(dh, jh, 'D-J')

    def _identify_v_j_junction(self):
        if 'dh' not in self.primary:
            vh = None if 'vh' not in self.primary else self.primary['vh']
            jh = None if 'jh' not in self.primary else self.primary['jh']
            self._identify_junction(vh, jh, 'V-J')

    def _check_for_stop_codon(self):
        if self.effective.codon and self.framed:
            if '*' in self.effective.codon:
                self.head['premature'] = True

    def _identify_cdr3(self):
        def locate_cdr3_start(vh):
            position = None
            gene = self.pipeline.resolver.gene_fetch(vh.subject_id)
            if 'cdr3' in gene.body:
                gene_cdr3 = gene.body['cdr3']
                position = vh.query_start + vh.reference_end - gene_cdr3['reference end']
                
                # because the gene annotation considers the cdr3 to start after the Cycteine
                position -= 3
                
                if position < 0 or position > self.sequence.length - 3:
                    self.make_comment('cdr3 start position out of range {}'.format(position))
                    position = None
                else:
                    if self.reverse_complemented:
                        acid = self.sequence.reversed.nucleotide[position:position+3]
                    else:
                        acid = self.sequence.nucleotide[position:position+3]

                    if acid in self.configuration.enumeration['nucleic to amino'] and self.configuration.enumeration['nucleic to amino'][acid] == 'C':
                        self.head['conserved cdr3 c'] = True
                    else:
                        self.make_comment('framing cycteine not conserved {}'.format(acid))
            return position

        def locate_cdr3_end(jh):
            position = None
            gene = self.pipeline.resolver.gene_fetch(jh.subject_id)
            if 'cdr3' in gene.body:
                gene_cdr3 = gene.body['cdr3']
                position = jh.query_start + jh.reference_end - gene_cdr3['reference start']

                # because the gene annotation considers the cdr3 to end before the Tryptophan
                position += 3

                if position < 0 or position > self.sequence.length - 3:
                    self.make_comment('cdr3 end position out of range {}'.format(position))
                    position = None
                else:
                    if self.reverse_complemented:
                        acid = self.sequence.reversed.nucleotide[position-3:position]
                    else:
                        acid = self.sequence.nucleotide[position-3:position]

                    if acid in self.configuration.enumeration['nucleic to amino'] and self.configuration.enumeration['nucleic to amino'][acid] == 'W':
                        self.head['conserved cdr3 w'] = True
                    else:
                        self.make_comment('framing tryptophan not conserved {}'.format(acid))
            return position

        jh = None if 'jh' not in self.primary else self.primary['jh']
        vh = None if 'vh' not in self.primary else self.primary['vh']
        if vh is not None and jh is not None:
            start = locate_cdr3_start(vh)
            end = locate_cdr3_end(jh)
            if start is not None and end is not None and start < end:
                cdr3 = Hit(self, { 
                    'FLAG': self.FLAG,
                    'region': 'cdr3',
                    'subject id': 'cdr3',
                    'query start': start,
                    'query end': end
                })
                if cdr3.valid:
                    print(to_json_document(cdr3.query))
                    # print(to_json_document(cdr3))
                    self._add_hit(cdr3)
                    self._pick_hit(cdr3)
                    cdr3.charge
                    cdr3.weight

    def _identify_chewback(self):
        for hit in self.hit:
            hit.chew_3_prime
            hit.chew_5_prime

    def _check_for_productivity(self):
        jh = None if 'jh' not in self.primary else self.primary['jh']
        vh = None if 'vh' not in self.primary else self.primary['vh']
        cdr3 = None if 'cdr3' not in self.primary else self.primary['cdr3']
        if (vh is not None and jh is not None and \
            cdr3 is not None and self.conserved_cdr3_cycteine and self.conserved_cdr3_tryptophan and \
            self.in_frame and not self.premature):
            if vh.functionality == 'F':
                if jh.functionality == 'F':
                    self.head['productive'] = True
                else:
                    self.make_comment('leading jh {} is non functional {}'.format(jh.gene, jh.functionality))
            else:
                self.make_comment('leading vh {} is non functional {}'.format(vh.gene, vh.functionality))

    def _analyze_quality(self):
        for name, hit in self.primary.items():
            self.head['{} length'.format(name)] = hit.query.length
            hit.node['average phread'] = float(sum(hit.query.phred)) / float((hit.query.length))

        self.head['average phred'] = float(sum(self.effective.phred)) / float((self.effective.length))

    def invalidate(self, message):
        self.head['valid'] = False
        self.make_comment(message)

    def make_comment(self, message):
        if 'comment' not in self.body:
            self.body['comment'] = []

        self.body['comment'].append(message)
        self.log.info('sample %s : %s', self.id, message)

    def reverse(self):
        self.body['sequence'] = self.body['sequence'].reversed

    def view(self, request):
        diagram = Diagram(self.pipeline, self, request)
        print(diagram.draw())

    def info(self):
        print(to_json_document(self))

    @property
    def abundance(self):
        return self.head['abundance']

    @property
    def has_dh(self):
        return 'dh' in self.primary

    @property
    def observable(self):
        return 'cdr3' in self.primary and \
            self.effective.codon and \
            self.cdr3_sequence['query'].codon and \
            len(self.cdr3_sequence['query'].codon) > 2

    @property
    def observation(self):
        if self._observation is None:
            if self.observable:
                self._observation = []
                breakdown = {
                    'feature': {
                        'id': self.id,
                        'cdr3 length': float(self.primary['cdr3']['query'].length),
                        'cdr3 charge': float(self.primary['cdr3']['charge']),
                        'cdr3 weight': float(self.primary['cdr3']['weight']),
                    }
                }

                breakdown['feature']['abundance'] = self.abundance
                breakdown['feature']['cdr3 nucleotide'] = self.cdr3_sequence['query'].nucleotide
                breakdown['feature']['cdr3 codon'] = self.cdr3_sequence['query'].codon
                breakdown['feature']['effective nucleotide'] = self.effective.nucleotide
                breakdown['feature']['effective codon'] = self.effective.codon

                phred = sorted(self.effective.phred)
                breakdown['feature']['effective min phred'] = float(phred[0])
                breakdown['feature']['effective low 4 phred'] = numpy.mean(phred[0:4])

                phred = sorted(self.cdr3_sequence['query'].phred)
                breakdown['feature']['cdr3 min phred'] = float(phred[0])
                breakdown['feature']['cdr3 low 4 phred'] = numpy.mean(phred[0:4])

                if self.has_dh:
                    breakdown['region'] = { 'vh': [], 'dh': [], 'jh': [] }
                    for hit in self.hit:
                        if hit.picked and hit.region in breakdown['region'].keys():
                            breakdown['region'][hit.region].append(hit)
                    breakdown['combination'] = len(breakdown['region']['vh']) * len(breakdown['region']['dh']) * len(breakdown['region']['jh'])
                    breakdown['feature']['weight'] = 1.0 / float(breakdown['combination'])

                    for feature in [
                        'chew',
                        'v3 chew',
                        'd5 chew',
                        'd3 chew',
                        'j5 chew',
                        'cdr3 length',
                        'cdr3 charge',
                        'cdr3 weight',
                        'vd length',
                        'vd n count',
                        'vd p count',
                        'dj length',
                        'dj n count',
                        'dj p count',
                        'n count',
                        'p count',
                    ]:
                        if feature not in breakdown['feature']:
                            breakdown['feature'][feature] = 0.0

                    # chewback
                    chew = [ v.chew_3_prime.length for v in breakdown['region']['vh'] if v.chew_3_prime is not None ]
                    if chew: breakdown['feature']['v3 chew'] = float(numpy.mean(chew))
                    chew = [ j.chew_5_prime.length for j in breakdown['region']['jh'] if j.chew_5_prime is not None ]
                    if chew: breakdown['feature']['j5 chew'] = float(numpy.mean(chew))
                    chew = [ d.chew_3_prime.length for d in breakdown['region']['dh'] if d.chew_3_prime is not None ]
                    if chew: breakdown['feature']['d3 chew'] = float(numpy.mean(chew))
                    chew = [ d.chew_5_prime.length for d in breakdown['region']['dh'] if d.chew_5_prime is not None ]
                    if chew: breakdown['feature']['d5 chew'] = float(numpy.mean(chew))
                    breakdown['feature']['chew'] = \
                        breakdown['feature']['v3 chew'] + \
                        breakdown['feature']['j5 chew'] + \
                        breakdown['feature']['d3 chew'] + \
                        breakdown['feature']['d5 chew']

                    # junction
                    if 'V-D' in self.primary:
                        breakdown['feature']['vd length'] = float(self.primary['V-D'].query.length)
                        breakdown['feature']['vd n count'] = float(self.primary['V-D'].palindrome.count('N'))
                        breakdown['feature']['vd p count'] = float(self.primary['V-D'].palindrome.count('P'))
                        breakdown['feature']['n count'] += breakdown['feature']['vd n count']
                        breakdown['feature']['p count'] += breakdown['feature']['vd p count']
                    if 'D-J' in self.primary:
                        breakdown['feature']['dj length'] = float(self.primary['D-J'].query.length)
                        breakdown['feature']['dj n count'] = float(self.primary['D-J'].palindrome.count('N'))
                        breakdown['feature']['dj p count'] = float(self.primary['D-J'].palindrome.count('P'))
                        breakdown['feature']['n count'] += breakdown['feature']['dj n count']
                        breakdown['feature']['p count'] += breakdown['feature']['dj p count']

                    for vh in breakdown['region']['vh']:
                        for dh in breakdown['region']['dh']:
                            for jh in breakdown['region']['jh']:
                                observation = deepcopy(breakdown['feature'])
                                observation['vh'] = vh.gene
                                observation['dh'] = dh.gene
                                observation['jh'] = jh.gene
                                observation['vh family'] = vh.family
                                observation['dh family'] = dh.family
                                self._observation.append(observation)
                else:
                    breakdown['region'] = { 'vh': [], 'jh': [] }
                    for hit in self.hit:
                        if hit.picked and hit.region in breakdown['region'].keys():
                            breakdown['region'][hit.region].append(hit)
                    breakdown['combination'] = len(breakdown['region']['vh']) * len(breakdown['region']['jh'])
                    breakdown['feature']['weight'] = 1.0 / float(breakdown['combination'])

                    for feature in [
                        'chew',
                        'v3 chew',
                        'j5 chew',
                        'cdr3 length',
                        'cdr3 charge',
                        'cdr3 weight',
                        'vj length',
                        'vj n count',
                        'vj p count',
                        'n count',
                        'p count',
                     ]:
                        if feature not in breakdown['feature']:
                            breakdown['feature'][feature] = 0.0

                    # chewback
                    chew = [ float(j.chew_5_prime.length) for j in breakdown['region']['jh'] if j.chew_5_prime is not None ]
                    if chew: breakdown['feature']['j5 chew'] = float(numpy.mean(chew))
                    chew = [ float(v.chew_3_prime.length) for v in breakdown['region']['vh'] if v.chew_3_prime is not None ]
                    if chew: breakdown['feature']['v3 chew'] = float(numpy.mean(chew))
                    breakdown['feature']['chew'] = \
                        breakdown['feature']['j5 chew'] + \
                        breakdown['feature']['v3 chew']

                    # junction
                    if 'V-J' in self.primary:
                        breakdown['feature']['vj length'] = float(self.primary['V-J'].query.length)
                        breakdown['feature']['vj n count'] = float(self.primary['V-J'].palindrome.count('N'))
                        breakdown['feature']['vj p count'] = float(self.primary['V-J'].palindrome.count('P'))

                    for vh in breakdown['region']['vh']:
                        for jh in breakdown['region']['jh']:
                            observation = deepcopy(breakdown['feature'])
                            observation['vh'] = vh.gene
                            observation['jh'] = jh.gene
                            observation['vh family'] = vh.family
                            self._observation.append(observation)

        return self._observation

    @property
    def cdr3_sequence(self):
        if 'cdr3' in self.primary:
            return self.primary['cdr3']
        else:
            return None

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def head(self):
        if 'head' not in self.node:
            self.node['head'] = {}
        return self.node['head']

    @property
    def body(self):
        if 'body' not in self.node:
            self.node['body'] = {}
        return self.node['body']

    @property
    def key(self):
        if 'id sha1' not in self.head and self.id is not None:
            self.head['id sha1'] = hashlib.sha1(self.id.encode('utf8')).hexdigest()
        return self.head['id sha1']

    @property
    def document(self):
        return to_document(self.node)

    @property
    def fasta(self):
        return to_fasta(self.id, self.sequence.nucleotide, None, self.configuration.constant['fasta line length'])

    @property
    def fastq(self):
        return to_fastq(self.id, self.sequence.nucleotide, self.sequence.quality)

    @property
    def id(self):
        return self.head['id']

    @property
    def length(self):
        return self.sequence.length

    @property
    def FLAG(self):
        return self.head['FLAG']

    @FLAG.setter
    def FLAG(self, value):
        self.head['FLAG'] = value

    @property
    def reverse_complemented(self):
        return self.FLAG & 0x10 == 0x10

    @reverse_complemented.setter
    def reverse_complemented(self, value):
        if value:
            self.FLAG |= 0x10
        else:
            self.FLAG &= ~0x10

    @property
    def library(self):
        return self.head['library']

    @property
    def sequence(self):
        return self.body['sequence']

    @property
    def effective(self):
        if 'effective' not in self.body:
            # find the borders
            start = self.sequence.length
            end = 0
            for hit in self.hit:
                if hit.valid and hit.type == 'gene segment':
                    if self.reverse_complemented:
                        start = min(start, self.sequence.length - hit.query_end)
                        end = max(end, self.sequence.length - hit.query_start)
                    else:
                        start = min(start, hit.query_start)
                        end = max(end, hit.query_end)

            if start > end:
                # no hits
                if self.reverse_complemented:
                    self.body['effective'] = self.sequence.reversed
                else:
                    self.body['effective'] = self.sequence
                self.body['effective query start'] = 0
                self.body['effective query end'] = self.sequence.length
            else:
                self.body['effective query start'] = start
                self.body['effective query end'] = end
                if self.reverse_complemented:
                    self.body['effective'] = self.sequence.reversed.crop(start, end)
                else:
                    self.body['effective'] = self.sequence.crop(start, end)
            # print(to_json_document(self.body['effective']))
        return self.body['effective']

    @property
    def valid(self):
        return self.head['valid']

    @property
    def hit(self):
        if 'hit' not in self.body:
            self.body['hit'] = []
        return self.body['hit']

    @property
    def primary(self):
        if self._primary is None:
            self._primary = {}
            if 'primary' in self.body:
                for k,v in self.body['primary'].items():
                    self._primary[k] = self.index[v]
        return self._primary

    @property
    def comment(self):
        if 'comment' in self.body:
            return self.body['comment']
        else:
            return None

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
    def conserved_cdr3_cycteine(self):
        return self.head['conserved cdr3 c']

    @property
    def conserved_cdr3_tryptophan(self):
        return self.head['conserved cdr3 w']

    @property
    def productive(self):
        return self.head['productive']

class Pivot(object):
    def __init__(self, study, node=None, name=None, rotate=None):
        self.log = logging.getLogger('Pivot')
        self.study = study
        self.node = node
        if self.node is None:
            self.node = {
                'name': name,
                'axis': None,
                'rotate': None,
                'pivot': None,
                'count': 0,
                'weight': [],
                'residual': {}
            }
            for feature in self.row:
                # TODO types of features
                self.residual[feature['key']] = dict([(column['key'], None) for column in feature['column']])
                self.residual[feature['key']]['value'] = []

            if rotate:
                self.node['axis'] = rotate[0]
                self.node['pivot'] = {}
                if len(rotate) > 1:
                    self.node['rotate'] = rotate[1:]
        else:
            if isinstance(self.node['pivot'], list):
                self.node['pivot'] = dict([(node['name'], Pivot(self.study, node)) for node in self.node['pivot']])

    @property
    def configuration(self):
        return self.study.configuration

    @property
    def empty(self):
        return not self.count > 0

    @property
    def document(self):
        document = to_document(self.node)
        if isinstance(document['pivot'], dict):
            document['pivot'] = list(document['pivot'].values())
        return document

    @property
    def row(self):
        return self.study.row

    @property
    def feature(self):
        return self.study.feature

    @property
    def axis(self):
        return self.node['axis']

    @property
    def name(self):
        return self.node['name']

    @property
    def rotate(self):
        return self.node['rotate']

    @property
    def pivot(self):
        return self.node['pivot']

    @property
    def empty(self):
        return not self.count > 0

    @property
    def count(self):
        return self.node['count']

    @count.setter
    def count(self, value):
        self.node['count'] = value

    @property
    def weight(self):
        return self.node['weight']

    @property
    def residual(self):
        return self.node['residual']

    def add_observation(self, observation):
        self.count += observation['weight']
        self.weight.append(observation['weight'] * observation['abundance'])
        for key,value in self.residual.items():
            value['value'].append(observation[key])

        if self.axis is not None:
            name = observation[self.axis]
            if name not in self.pivot:
                self.pivot[name] = Pivot(self.study, None, name, self.rotate)
            pivot = self.pivot[name]
            pivot.add_observation(observation)

    def finalize(self):
        if not self.empty:
            self.node['weight'] = numpy.array(self.node['weight'])
            for key,feature in self.residual.items():
                feature['value'] = numpy.array(feature['value'])
                feature['min'] = numpy.amin(feature['value'])
                feature['max'] = numpy.amax(feature['value'])
                feature['mean'] = numpy.average(feature['value'], weights=self.weight)
                feature['std'] = math.sqrt(numpy.average(pow(feature['value'] - feature['mean'], 2), weights=self.weight))

        if self.axis is not None:
            for pivot in self.pivot.values():
                pivot.finalize()

    def collect(self, template, buffer):
        if 'depth' not in template:
            template['depth'] = 0

        if 'rotate' not in template:
            template['rotate'] = []

        if template['depth'] > 0:
            for name, pivot in self.pivot.items():
                iteration = {
                    'depth': template['depth'] - 1,
                    'rotate': template['rotate'] + [ name ]
                }
                pivot.collect(iteration, buffer)
        else:
            template['count'] = round(self.count, self.configuration.constant['float accuracy'])
            template['weight'] = round(numpy.sum(self.weight), self.configuration.constant['float accuracy'])
            template['residual'] = self.residual
            buffer.append(template)
        return buffer

class Study(object):
    def __init__(self, pipeline, node=None, request=None):
        self.log = logging.getLogger('Study')
        self.pipeline = pipeline
        self.node = node
        self.load(request)

    def load(self, request=None):
        if self.node is not None:
            if 'body' not in self.node:
                self.node['body'] = {}

            if 'pivot' in self.body:
                self.body['pivot'] = Pivot(self, self.body['pivot'], None, self.rotate)

        elif request is not None:
            self.node = {
                'head' :{
                    'id': request.hash_for('study'),
                    'request': request.document_for('study')
                },
                'body': { }
            }
            self.head['name'] = request.name if request.name else self.id
            self.reset()
        else:
            raise ValueError('insufficient information to create study')

    def delete(self):
        path = os.path.join(self.configuration.cache['path'], self.id)
        if os.path.exists(path):
            os.remove(path)

    def persist(self):
        path = os.path.join(self.configuration.cache['path'], self.id)
        try:
            with io.open(path, 'wb') as file:
                pickle.dump(self.document, file)
            self.log.debug('pickled study %s to %s', self.id, path)
        except OSError as e:
            self.log.error('failed writing pickled study %s to %s', self.id, path)
            self.log.debug(str(e))

    def restore(self):
        path = os.path.join(self.configuration.cache['path'], self.id)
        if os.path.exists(path):
            with io.open(path, 'rb') as file:
                self.node = pickle.load(file)
                self.load()
        else:
            raise ValueError('missing cached object for {}'.format(self.id))

    def reset(self):
        self.count = 0
        self.body['repertoire'] = dict([(pivot, {}) for pivot in self.rotate])
        if self.request['repertoire']:
            for pivot in self.request['repertoire'].values():
                for predicate in pivot:
                    if predicate['pivot'] in self.repertoire:
                        self.repertoire[predicate['pivot']][predicate['key']] = deepcopy(predicate)
                        self.repertoire[predicate['pivot']][predicate['key']]['weight'] = 0.0

    def __str__(self):
        return '\t'.join([self.id, self.name, str(self.count)])

    @property
    def document(self):
        return to_document(self.node)

    @property
    def persisted(self):
        return '_id' in self.node and self.node['_id']

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def head(self):
        return self.node['head']

    @property
    def body(self):
        return self.node['body']

    @property
    def id(self):
        return self.head['id']

    @property
    def name(self):
        return self.head['name']

    @property
    def request(self):
        return self.head['request']

    @property
    def preset(self):
        return self.request['preset']

    @property
    def unique(self):
        return self.preset['unique']

    @property
    def row(self):
        return self.preset['row']

    @property
    def feature(self):
        return self.preset['feature']

    @property
    def rotate(self):
        return self.preset['rotate']

    @property
    def pivot(self):
        if 'pivot' not in self.body:
            if self.persisted:
                self.restore()
            else:
                self.body['pivot'] = Pivot(self, None, None, self.rotate)

        return self.body['pivot']

    @property
    def repertoire(self):
        if 'repertoire' not in self.body:
            if self.persisted:
                self.restore()
        return self.body['repertoire']

    @property
    def empty(self):
        return not self.count > 0

    @property
    def count(self):
        return self.head['count']

    @count.setter
    def count(self, value):
        self.head['count'] = value

    def index_sample(self, sample):
        self.count += 1
        for pivot in self.rotate:
            for observation in sample.observation:
                if observation[pivot] not in self.repertoire[pivot]:
                    self.repertoire[pivot][observation[pivot]] = {
                        'key': observation[pivot],
                        'name': observation[pivot],
                        'weight': 0.0,
                    }
                self.repertoire[pivot][observation[pivot]]['weight'] += observation['weight']

    def add_sample(self, sample):
        taken = True
        if sample.observable:
            if not self.unique or len(sample.observation) == 1:
                # check the observations
                for observation in sample.observation:
                    for pivot in self.rotate:
                        if pivot not in observation:
                            self.log.debug('ignoring sample %s with missing pivot %s', sample.id, pivot)
                            taken = False
                            break

                    if not taken:
                        break

                    for feature in self.row:
                        if feature['key'] not in observation:
                            self.log.debug('ignoring sample %s with missing feature %s', sample.id, feature['key'])
                            taken = False
                            break

                    if not taken:
                        break

                if taken:
                    self.index_sample(sample)
                    for observation in sample.observation:
                        self.pivot.add_observation(observation)
            else:
                self.log.debug('ignoring sample %s with multiple observations', sample.id)
        else:
            self.log.debug('ignoring sample %s with no observations', sample.id)

    def finalize(self):
        self.pivot.finalize()

    def csv(self, request):
        if self.pivot is not None:
            buffer = []
            preset = { 
                'depth': request.instruction['depth'],
            }
            if request.preset:
                preset['row'] = request.preset['row']
            print(to_csv(self.encode_csv_head(preset)))

            self.encode_csv(preset, buffer)
            for row in buffer:
                print(to_csv(row))

    def json(self, request):
        self.pivot
        print(to_json_document(self))

    def collect(self, template, buffer):
        return self.pivot.collect({ 'depth': template['depth'] }, buffer)

    def encode_csv_head(self, preset=None):
        # if preset is not None and 'row' in preset:
        #     template = { 'row': deepcopy(preset['row']) }
        # else:
        #     template = { 'row': deepcopy(self.row) }

        template = { 'row': deepcopy(self.row) }
        if preset:
            template.update(preset)

        if 'depth' not in template or template['depth'] is None:
            template['depth'] = len(self.rotate)

        template['depth'] = min(template['depth'], len(self.rotate))

        head = self.rotate[:template['depth']] + [ 'study', 'weight', 'count']
        for feature in template['row']:
            for column in feature['column']:
                head.append('{}_{}'.format(feature['key'].replace(' ', '_'), column['key'].replace(' ', '_')))
        return head

    def encode_csv(self, preset, buffer):
        def construct(template):
            row = self.rotate[:template['depth']] + [ 'study', 'weight', 'count' ]
            for feature in template['row']:
                for column in feature['column']:
                    row.append('{} {}'.format(feature['key'], column['key']))
            structure = {
                'head': row,
                'index': {},
                'name': {},
            }
            for index,name in enumerate(row):
                structure['index'][index] = name
                structure['name'][name] = index
            return structure

        def table(template, collection, buffer):
            weight = float(self.count)
            for record in collection:
                row = record['rotate'] + [ self.name, record['weight'], record['count'] ]
                for feature in template['row']:
                    for column in feature['column']:
                        if feature['key'] in record['residual'] and column['key'] in record['residual'][feature['key']]:
                            row.append(round(record['residual'][feature['key']][column['key']], self.configuration.constant['float accuracy']))
                        else:
                            row.append(None)
                buffer.append(row)

        def sort(template, structure, buffer):
            for sort in template['sort']:
                if sort['column'] in structure['name']:
                    index = structure['name'][sort['column']]
                buffer.sort(key=lambda x: x[index], reverse=sort['reverse'])

        template = {
            'row': deepcopy(self.row),
            'sort': deepcopy(self.preset['sort'])
        }

        if preset:
            template.update(preset)

        if 'depth' not in template or template['depth'] is None:
            template['depth'] = len(self.rotate)
        template['depth'] = min(template['depth'], len(self.rotate))

        structure = construct(template)
        collection = self.collect(template, [])
        table(template, collection, buffer)
        sort(template, structure, buffer)

        return buffer

    def tiled_pivot_expression(self, name, interval):
        pivot = sorted(list(self.repertoire[name].values()), key=lambda x: x['start'])
        intensity = {
            'max name length': max([len(predicate['name']) for predicate in pivot]),
            'pivot count': len(pivot),
            'tile length': interval,
            'start': pivot[0]['start'],
            'end': pivot[-1]['end'],
            'tiles': [],
        }
        intensity['length'] = intensity['end'] - intensity['start']
        intensity['tile count'] = math.ceil(float(intensity['length']) / float(intensity['tile length']))
        intensity['tiled length'] = intensity['tile count'] * intensity['tile length']
        intensity['tiled start'] = intensity['start'] - math.floor(float(intensity['tiled length'] - intensity['length']) / 2.0)
        intensity['tiled end'] = intensity['tiled start'] + intensity['tiled length']

        for index in range(intensity['tile count']):
            intensity['tiles'].append(
                {
                    'intensity': 0.0,
                    'pivot':[],
                    'start': intensity['tiled start'] + intensity['tile length'] * index,
                    'end': intensity['tiled end'] + intensity['tile length'] * (index + 1),
                }
            )

        for predicate in pivot:
            tile = None
            left = math.floor(float(predicate['start'] - intensity['tiled start']) / float(intensity['tile length']))
            right = math.floor(float(predicate['end'] - intensity['tiled start']) / float(intensity['tile length']))
            if left == right:
                tile = intensity['tiles'][left]
            else:
                split = intensity['tiled start'] + ( right * intensity['tile length'] )
                left_part = float(split - predicate['start']) / float(predicate['length'])
                right_part = float(predicate['end'] - split) / float(predicate['length'])
                if left_part > right_part:
                    tile = intensity['tiles'][left]
                else:
                    tile = intensity['tiles'][right]
            # tile['intensity'] += 1
            tile['intensity'] += predicate['weight']
            tile['pivot'].append(predicate)
        return intensity

    def tiled_expression(self, buffer):
        density = {
            'vh': {
                'order': 3,
                'interval': 70000,
            },
            'dh': {
                'order': 2,
                'interval': 2500,
            },
            'jh': {
                'order': 1,
                'interval': 400,
            }
        }
        for pivot, node in density.items():
            if pivot in self.rotate:
                node['intensity'] = self.tiled_pivot_expression(pivot, node['interval'])
                node['total'] = float(sum([tile['intensity'] for tile in node['intensity']['tiles']]))
                node['factor'] = float(node['intensity']['pivot count']) / node['total']
                for tile in node['intensity']['tiles']:
                    row = {
                        'study': self.name,
                        'pivot': pivot,
                        'start': tile['start'],
                        'end': tile['end'],
                        'intensity': tile['intensity'],
                        'adjusted': float(tile['intensity'] * node['factor']),
                        'element': ';'.join([e['name'] for e in tile['pivot']])
                    }
                    buffer.append(row)

    def info(self):
        print(to_json_document(self))

class Diagram(object):
    def __init__(self, pipeline, sample, request):
        self.log = logging.getLogger('Diagram')
        self.pipeline = pipeline
        self.sample = sample
        self.sequence = self.sample.effective
        self.offset = self.sample.body['effective query start']

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
        profile = request.profile_for('sample')
        if 'format' in profile and 'diagram' in profile['format']:
            if 'track' in profile['format']['diagram']:
                for k,v in profile['format']['diagram']['track'].items():
                    self.query[k] = v

        for track in sample.hit:
            if track['region'] not in self.configuration.region.keys() or all([k in track and track[k] == v for k,v in self.query.items()]):
                self.track.append(track)
                
        for k in profile['format']['diagram']['feature']:
            if k in self.configuration.diagram['prototype']:
                feature = dict(self.configuration.diagram['prototype'][k])
                if feature['width'] == 'auto':
                    feature['width'] = 0
                    for track in self.track:
                        if feature['value'] in track:
                            feature['width'] = max(feature['width'], len(feature['format'](track[feature['value']])))
                self.pattern['feature'].append(feature)
                self.pattern['phrase'].append('{{: <{:}}}'.format(feature['width']))
                self.pattern['title'].append(feature['title'])
                self.pattern['width'] += max(feature['width'], len(feature['title']))
        self.pattern['width'] += (len(self.pattern['phrase']) - 1)
        self.pattern['phrase'] = ' '.join(self.pattern['phrase'])
        self.width['sample read frame'] = self.sequence.read_frame
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

    def _draw_comment(self, buffer):
        if self.sample.comment:
            buffer.write('{: <{}}'.format('Comment', self.width['diagram start']))
            buffer.write(' | '.join(self.sample.comment))
            buffer.write('\n')

    def _draw_summary(self, buffer):
        b = []
        b.append(self.sample.library)
        b.append(str(self.sample.abundance))
        b.append(self.sample.id)
        if 'vh' in self.sample.primary:
            vh = self.sample.primary['vh']
            b.append('{} : {}'.format(vh.region, vh.gene))

        if 'dh' in self.sample.primary:
            dh = self.sample.primary['dh']
            b.append('{} : {}'.format(dh.region, dh.gene))

        if 'jh' in self.sample.primary:
            jh = self.sample.primary['jh']
            b.append('{} : {}'.format(jh.region, jh.gene))

        if 'cdr3' in self.sample.primary:
            cdr3 = self.sample.primary['cdr3']
            if 'charge' in cdr3 and 'weight' in cdr3:
                b.append('{} : {:.3} {}'.format(cdr3.region, cdr3.node['charge'], cdr3.node['weight']))

        if self.sample.productive:
            b.append('productive')
        else:
            if self.sample.framed:
                b.append('in frame' if self.sample.in_frame else 'out of frame')
                
            if self.sample.premature:
                b.append('premature')

        if 'average phred' in self.sample.head:
            b.append('Q {:.4}'.format(self.sample.head['average phred']))
        
        buffer.write(' | '.join(b))
        buffer.write('\n')

    def _draw_track_title(self, buffer):
        buffer.write(self.pattern['phrase'].format(*(self.pattern['title'])))

    def _draw_coordinates(self, buffer):
        buffer.write(' ' * self.width['diagram start'])
        if self.width['sample read frame'] > 0:
            buffer.write('{: <{}}'.format(self.offset, self.width['sample read frame']))
            buffer.write(self.gap)
            
        for index in range(self.width['sample read frame'] + self.offset, self.sequence.length + self.offset, 3):
            buffer.write('{: <3}'.format(index))
            buffer.write(self.gap)
        buffer.write('\n')

    def _draw_nucleotide_sequence(self, buffer):
        buffer.write(' ' * self.width['diagram start'])
        if self.width['sample read frame'] > 0:
            buffer.write(self.sequence.nucleotide[0:self.width['sample read frame']])
            buffer.write(self.gap)

        for index in range(self.width['sample read frame'], self.sequence.length, 3):
            buffer.write(self.sequence.nucleotide[index:index + 3])
            buffer.write(self.gap)
        buffer.write('\n')

    def _draw_quality_sequence(self, buffer):
        buffer.write(' ' * self.width['diagram start'])
        if self.width['sample read frame'] > 0:
            buffer.write(self.sequence.quality[0:self.width['sample read frame']])
            buffer.write(self.gap)

        for index in range(self.width['sample read frame'], self.sequence.length, 3):
            buffer.write(self.sequence.quality[index:index + 3])
            buffer.write(self.gap)
        buffer.write('\n')

    def _draw_codon_sequence(self, buffer):
        if self.sample.framed:
            buffer.write(' ' * self.width['diagram start'])
            if self.width['sample read frame'] > 0:
                buffer.write(' ' * self.width['sample read frame'])
                buffer.write(self.gap)

            for codon in self.sequence.codon:
                buffer.write('{: <3}'.format(codon))
                buffer.write(self.gap)
            buffer.write('\n')

    def _draw_track_without_subject(self, buffer, track):
        offset = self.node['track offset'][track['uuid']]
        if 'subject' not in track and 'query' in track:
            if offset > 0: buffer.write(' ' * offset)
            if track['query'].read_frame > 0:
                buffer.write('-' * min(track['query'].read_frame, track['query'].length))
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
                    buffer.write(track['palindrome'][0:min(track['query'].read_frame, track['query'].length)])
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
        start = track['query start'] - self.offset
        if start > 0:
            result = start
            if self.sample.framed:
                if self.width['sample read frame'] > 0 and start >= self.width['sample read frame']:
                    result += int((start - self.width['sample read frame']) / 3) + 1
                else:
                    result += int(start / 3)
        return result

    def _draw_track(self, buffer, track):
        buffer.write(self.format_track_features(track))
        buffer.write(' ' * self.width['table padding'])
        self._draw_track_without_subject(buffer, track)
        self._draw_track_with_subject(buffer, track)

    def draw(self):
        buffer = StringIO()
        buffer.write('\n')
        self._draw_comment(buffer)
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
    def __init__(self, pipeline, request):
        self.log = logging.getLogger('Block')
        self.pipeline = pipeline
        self.node = {
            'capacity': self.configuration.constant['buffer size'],
            'preset': request.preset_for('sample'),
            'library': request.instruction['library'],
            'buffer': None,
            'lookup': None
        }
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
    def capacity(self):
        return self.node['capacity']

    @property
    def size(self):
        return len(self.buffer)

    @property
    def empty(self):
        return self.size == 0

    @property
    def library(self):
        return self.node['library']

    @property
    def preset(self):
        return self.node['preset']

    @property
    def organism(self):
        return self.configuration.selected['organism']

    @property
    def strain(self):
        return self.configuration.selected['strain']

    @property
    def chain(self):
        return self.configuration.selected['chain']

    @property
    def region(self):
        return self.configuration.selected['chain']['region']

    @property
    def full(self):
        return not self.size < self.capacity

    def reset(self):
        self.node['buffer'] = []
        self.node['lookup'] = {}

    def add(self, sample):
        #if sample is not None and sample.sequence.length > 0:
        if sample.id not in self.lookup:
            self.buffer.append(sample)
            self.lookup[sample.id] = sample
        else:
            self.log.error('%s already present', sample.id)

    def fill(self):
        if self.read():
            self.search()
            for sample in self.buffer:
                sample.analyze()
                # print(to_json_document(sample))
                # print(sample.productive)
        return not self.empty

    def read(self):
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
                            sample = Sample(self.pipeline, id=line[1:], library=self.library)
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
        for key in self.chain['alignment order']:
            self.bwa(self.region[key])
            self.determine_region(key)

    def determine_region(self, region):
        for sample in self.buffer:
            sample.determine_region(region)

    def bwa(self, region):
        command = Command(self.pipeline, 'bwa mem', region['bwa mem'])
        command.args.append(region['bwa fasta path'])
        command.args.append('-')
        output, error = command.process.communicate(input=self.to_fastq().read().encode('utf8'))
        if output:
            buffer = StringIO(output.decode('utf8'))
            state = 0
            for line in buffer:
                line = line.strip('\n')
                if line[0] != '@':
                    hit = Hit.from_bwa(self, line, region['key'], region['type'])
                    if hit and hit.valid:
                        hit.sample.add_bwa_hit(hit)

    def igblast(self):
        chain = self.configuration.selected['chain']
        Command = command(self.pipeline, 'igblast', chain['igblast'])
        command.cwd = chain['igblast database directory']
        output, error = command.process.communicate(input=self.to_fasta().read().encode('utf8'))
        if output:
            buffer = StringIO(output.decode('utf8'))
            buffer.seek(0)
            state = 0
            sample = None
            for line in buffer:
                line = line.strip()
                
                if line.startswith('# IGBLASTN '):
                    # this is the start of a new record
                    if state == 2:
                        sample.invalidate('no alignments found')
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
                        
                elif state == 2 and line.startswith(self.configuration.expression['igblast reversed query']):
                    # this means the hits will be for the reverse complement strand
                    sample.reverse()
                    
                elif state == 2 and line.startswith('# Hit table '):
                    # this means the hit table is on the next line
                    state = 3
                    
                elif state == 3 and not line.startswith('#'):
                    hit = Hit.from_igblast(self, line, sample.id)
                    if hit:
                        hit.sample.add_igblast_hit(hit)

    def to_fastq(self):
        buffer = StringIO()
        for sample in self.buffer:
            residual = sample.residual
            if residual is not None:
                buffer.write(to_fastq(sample.id, residual.nucleotide, residual.quality))
                buffer.write('\n')
        buffer.seek(0)
        return buffer

    def to_fasta(self, query=None):
        buffer = StringIO()
        for sample in self.buffer:
            if query is None or all([k in sample.head and sample.head[k] == v for k,v in query.items()]):
                buffer.write(to_fasta(sample.id, sample.sequence.nucleotide, None, self.configuration.constant['fasta line length']))
                buffer.write('\n')
        buffer.seek(0)
        return buffer

class Library(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Library')
        self.pipeline = pipeline
        self.node = node

    def __str__(self):
        return '| ' + ' | '.join(['{:<9}', '{:<3}', '{:<2}', '{:<2}', '{:<30}', '{:12}', '{:8}', '{:8}', '{:<25}']).format(
            '' if 'strain' not in self.head else str(self.head['strain']), 
            '' if 'exposure' not in self.head else str(self.head['exposure']), 
            '' if 'biological repetition' not in self.body else str(self.body['biological repetition']), 
            '' if 'technical repetition' not in self.body else str(self.body['technical repetition']), 
            '' if 'tissue' not in self.head else str(self.head['tissue']),
            str(self.count), 
            str(self.valid_count), 
            '0' if self.valid_count == 0 else '{:.4}'.format(self.productive_count / self.valid_count * 100.0), 
            self.id
        )

    @property
    def document(self):
        return to_document(self.node)

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def id(self):
        return self.head['id']

    @property
    def head(self):
        return self.node['head']

    @property
    def body(self):
        return self.node['body']

    @property
    def strain(self):
        return None if 'strain' not in self.head else self.head['strain']

    @property
    def reference(self):
        if 'reference' in self.body:
            return self.body['reference']
        else:
            return None

    @property
    def count(self):
        if 'count' in self.body:
            return self.body['count']
        else:
            return None

    @count.setter
    def count(self, value):
        if value is None and 'count' in self.body:
            del self.body['count']
        else:
            self.body['count'] = value

    @property
    def valid_count(self):
        if 'valid count' in self.body:
            return self.body['valid count']
        else:
            return None

    @valid_count.setter
    def valid_count(self, value):
        if value is None and 'valid count' in self.body:
            del self.body['valid count']
        else:
            self.body['valid count'] = value

    @property
    def productive_count(self):
        if 'count productive' in self.body:
            return self.body['count productive']
        else:
            return None

    @productive_count.setter
    def productive_count(self, value):
        if value is None and 'count productive' in self.body:
            del self.body['count productive']
        else:
            self.body['count productive'] = value

class FastqFilteringReader(object):
    def __init__(self, pipeline, request=None):
        self.log = logging.getLogger('Filter')
        self.pipeline = pipeline
        self.input = None
        self.output = None
        self.forward = None
        self.reverse = None
        self.merged = None
        self.node = {
            'eof': False,
            'limit': None,
            'trim start': False,
            'trim end': True,
            'forward path': None,
            'reverse path': None,
            'trimming window': 3,
            'trimming error threshold': 0.5,
            'adapter mismatch': 3,
            'adapter overlap': 6,
            'adapter score threshold': 40,
            'minimum read length': 130,
            'error threshold': 1,
            'merge overlap threshold': 70,
            'merge score threshold': 450,
            'quota': 32,
            'paired': False,
            'filtering': True,
            'reverse adapter scan': True,
            'adapter': [
                'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                'TACACTCTTTCCCTACACGACGCTCTTCCGATCT',
            ],
            'retain short': True,
            'retain unmerged': True,
            'substitution quality correction': 28,
            'statistic': {
                'total': 0,
                'matched adapter': 0,
                'adapter mismatch': {},
                'adapter overlap': {},
                'adapter score': {},
                'merged': 0,
                'merged mismatch': {},
                'merged overlap': {},
                'merged score': {},
                'merged offset': {},
            }
        }

        # load parameters from request
        if request is not None:
            for term in [
                'limit',
                'forward path',
                'reverse path',
            ]:
                if term in request.instruction:
                    self.node[term] = request.instruction[term]

        self.node['substitution quality correction'] = float(self.node['substitution quality correction']) * 2.0

        # make sure all reverse complements of adapters are present
        adapter = set()
        for sequence in self.node['adapter']:
            adapter.add(sequence)
            adapter.add(complement(sequence))
        self.node['adapter'] = list(adapter)

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def limit(self):
        return self.node['limit']

    @property
    def quota(self):
        return self.node['quota']

    @property
    def paired(self):
        return self.node['paired']

    @property
    def ambiguity(self):
        if 'ambiguity' not in self.node:
            self.node['ambiguity'] = {}
            for i in self.configuration.enumeration['iupac nucleic acid notation'].keys():
                for j in self.configuration.enumeration['iupac nucleic acid notation'].keys():
                    if i not in self.node['ambiguity']:
                        self.node['ambiguity'][i] = {}
                    self.node['ambiguity'][i][j] = set(self.configuration.enumeration['iupac nucleic acid notation'][i]['option']).union(set(self.configuration.enumeration['iupac nucleic acid notation'][j]['option']))
            for i in self.configuration.enumeration['iupac nucleic acid notation'].keys():
                for j in self.configuration.enumeration['iupac nucleic acid notation'].keys():
                    c = self.node['ambiguity'][i][j]
                    for k,v in self.configuration.enumeration['iupac nucleic acid notation'].items():
                        o = set(v['option'])
                        if o == c:
                            self.node['ambiguity'][i][j] = k
        return self.node['ambiguity']

    @property
    def filtering(self):
        return self.node['filtering']

    @property
    def statistic(self):
        return self.node['statistic']

    @property
    def count(self):
        return self.statistic['total']

    @property
    def size(self):
        return len(self.input)

    @property
    def enough(self):
        return not (self.limit is None or self.count < self.limit)

    @property
    def empty(self):
        return not self.size > 0

    @property
    def full(self):
        return self.size >= self.quota

    @property
    def adapter(self):
        return self.node['adapter']

    def open(self):
        self.input = []
        self.output = []
        if self.node['forward path'] is not None:
            self.forward = io.open(self.node['forward path'], 'r')
        else:
            self.forward = sys.stdin

        if self.node['reverse path'] is not None:
            self.reverse = io.open(self.node['reverse path'], 'r')

        if self.forward is not None and self.reverse is not None:
            self.node['paired'] = True

        self.merged = sys.stdout

    def close(self):
        if self.node['forward path'] is not None:
            self.forward.close()

        if self.node['reverse path'] is not None:
            self.reverse.close()

        self.input = None
        self.output = None

    def read(self, file):
        state = 0
        segment = None
        for line in file:
            if not line:
                break

            else:
                line = line.strip()
                if line:
                    if state == 0 and line[0] == '@':
                        segment = { 'id': line[1:] }
                        state = 1

                    elif state == 1:
                        segment['nucleotide'] = line
                        state = 2

                    elif state == 2:
                        state = 3

                    elif state == 3:
                        segment['quality'] = line
                        state = 4
                        break

        if  state == 4 and \
            segment is not None and \
            len(segment['nucleotide']) == len(segment['quality']):
                segment['length'] = len(segment['nucleotide'])
                segment['expected error'] = expected_error(segment['quality'])
                segment['effective start'] = 0
                segment['effective end'] = segment['length']
                segment['effective length'] = segment['length']
        else:
            segment = None
            self.node['eof'] = True

        return segment

    def write(self):
        self.merged.write('\n'.join([ to_fastq(s['id'], s['nucleotide'],s['quality']) for s in self.output ]))
        self.merged.write('\n')

    def fill(self):
        self.input = []
        self.output = []
        while not self.full and not self.node['eof']:
            sample = None
            forward = self.read(self.forward)
            if forward:
                forward['origin'] = 'forward'
                sample =  {
                    'forward': forward,
                    'id': forward['id'].split()[0],
                }

                if self.paired:
                    reverse = self.read(self.reverse)
                    if reverse:
                        reverse['origin'] = 'reverse'
                        sample['reverse'] = reverse
                        if forward['id'].split()[0] != reverse['id'].split()[0]:
                            raise ValueError('forward and reverse read don\'t match {} {}'.format(forward['id'], reverse['id']))
                    else:
                        raise ValueError('missing matching reverse read for {}'.format(forward['id']))

            if sample is not None:
                self.input.append(sample)

        if self.filtering:
            for sample in self.input:
                self.locate_adapter(sample['forward'])
                if self.paired:
                    self.locate_adapter(sample['reverse'])
                    self.merge_sample(sample)
                self.pick_sample(sample)
        else:
            for sample in self.input:
                self.place_in_output(sample['forward'])

        self.statistic['total'] += len(self.input)
        return not self.empty

    def locate_adapter(self, segment):
        taken = False
        for index, adapter in enumerate(self.adapter):
            if self.find_adapter(segment, adapter):
                self.promote_adapter(index)
                taken = True
                break

    def find_adapter(self, segment, adapter):
        taken = False
        if self.node['reverse adapter scan']:
            right = segment['length']
            left = right - len(adapter)
            while left >= 0:
                if self.match_adapter(segment, adapter, left, right):
                    taken = True
                    break
                right -= 1
                left -= 1
        else:
            left = 0
            right = len(adapter)
            while right <= len(segment['nucleotide']):
                if self.match_adapter(segment, adapter, left, right):
                    taken = True
                    break
                right += 1
                left += 1

        if not taken:
            taken = self.locate_partial_adapter(segment, adapter)

        return taken

    def promote_adapter(self, index):
        if index > 0 and index < len(self.adapter):
            adapter = self.adapter.pop(index)
            self.adapter.insert(0, adapter)

    def locate_partial_adapter(self, segment, adapter):
        taken = False
        end = len(segment['nucleotide'])
        right = min(end, len(adapter))
        left = end - right
        while end - left >= self.node['adapter overlap']:
            if self.match_adapter(segment, adapter[:right], left, end):
                taken = True
                break
            right -= 1
            left += 1
        return taken

    def match_adapter(self, segment, adapter, left, right):
        result = False
        if len(adapter) == (right - left):
            distance = 0
            score = 0
            length = len(adapter)
            sequence = segment['nucleotide'][left:right]
            quality = segment['quality'][left:right]
            if length >= self.node['adapter overlap']:
                for i in range(length):
                    score += self.score_base_alignment(adapter[i], sequence[i], quality[i], quality[i])
                    if adapter[i] != sequence[i]:
                        distance += 1
                        if distance > self.node['adapter mismatch']:
                            break

                if  distance <= self.node['adapter mismatch'] and \
                    score >= self.node['adapter score threshold']:
                    segment['adapter'] = {
                        'start': left,
                        'end': right,
                        'distance': distance,
                        'score': score,
                        'length': right - left,
                        'nucleotide': sequence,
                    }
                    self.report_adapter(segment, segment['adapter'])
                    segment['effective start'] = 0
                    segment['effective end'] = left
                    segment['effective length'] = left
                    result = True
        else:
            raise ValueError('compared sequences must be of equal length')

        return result

    def report_adapter(self, segment, adapter):
        self.statistic['matched adapter'] += 1

        if adapter['score'] not in self.statistic['adapter score']:
            self.statistic['adapter score'][adapter['score']] = 0
        self.statistic['adapter score'][adapter['score']] += 1

        if adapter['distance'] not in self.statistic['adapter mismatch']:
            self.statistic['adapter mismatch'][adapter['distance']] = 0
        self.statistic['adapter mismatch'][adapter['distance']] += 1

        if adapter['length'] not in self.statistic['adapter overlap']:
            self.statistic['adapter overlap'][adapter['length']] = 0
        self.statistic['adapter overlap'][adapter['length']] += 1

        self.log.debug('%s found adapter %s at position %s with %s mismatch', segment['id'], adapter['nucleotide'], adapter['start'], adapter['distance'])

    def report_merge(self, sample):
        self.statistic['merged'] += 1

        score = round(sample['merged']['score'] / 10) * 10
        if score not in self.statistic['merged score']:
            self.statistic['merged score'][score] = 0
        self.statistic['merged score'][score] += 1

        mismatch = round(sample['merged']['mismatch'] / 10) * 10
        if mismatch not in self.statistic['merged mismatch']:
            self.statistic['merged mismatch'][mismatch] = 0
        self.statistic['merged mismatch'][mismatch] += 1

        overlap = round(sample['merged']['overlap'] / 10) * 10
        if overlap not in self.statistic['merged overlap']:
            self.statistic['merged overlap'][overlap] = 0
        self.statistic['merged overlap'][overlap] += 1

        offset = round(sample['merged']['offset'] / 10) * 10
        if offset not in self.statistic['merged offset']:
            self.statistic['merged offset'][offset] = 0
        self.statistic['merged offset'][offset] += 1

        self.log.debug('%s merged at offset %s with score %s', sample['merged']['id'], sample['merged']['offset'], sample['merged']['score'])

    def merge_sample(self, sample):
        forward = sample['forward']
        reverse = self.complement_segment(sample['reverse'])
        left = self.align(forward, reverse)
        right = self.align(reverse, forward)
        picked = None
        if left is not None or right is not None:
            if left is None or right['score'] > left['score']:
                picked = right
                picked['direction'] = 'right'
            else:
                picked = left
                picked['direction'] = 'left'

        if picked is not None and picked['score'] >= self.node['merge score threshold']:
            top = picked['top']
            bottom = picked['bottom']
            merged = {
                'id': sample['id'],
                'nucleotide': [],
                'quality': [],
                'offset': picked['offset'],
                'score': picked['score'],
                'overlap':0,
                'mismatch': 0,
                'origin': 'merged',
            }
            if picked['direction'] == 'right':
                t = top['effective start'] + picked['offset']
                b = bottom['effective start']
            else:
                t = top['effective start']
                b = bottom['effective start'] - picked['offset']

            while t < top['effective end'] or b < bottom['effective end']:
                te = None
                be = None
                if t >= top['effective start'] and t < top['effective end']:
                    te = t

                if b >= bottom['effective start'] and b < bottom['effective end']:
                    be = b

                if te is None and be is not None:
                    merged['nucleotide'].append(bottom['nucleotide'][be])
                    merged['quality'].append(bottom['quality'][be])
                elif te is not None and be is None:
                    merged['nucleotide'].append(top['nucleotide'][te])
                    merged['quality'].append(top['quality'][te])
                elif te is not None and be is not None:
                    merged['overlap'] += 1
                    if top['nucleotide'][te] != bottom['nucleotide'][be]:
                        merged['mismatch'] += 1

                    if top['quality'][te] > bottom['quality'][be]:
                        merged['nucleotide'].append(top['nucleotide'][te])
                    elif top['quality'][te] < bottom['quality'][be]:
                        merged['nucleotide'].append(bottom['nucleotide'][be])
                    else:
                        nucleotide = self.ambiguate(top['nucleotide'][te], bottom['nucleotide'][be])
                        merged['nucleotide'].append(nucleotide)

                    merged['quality'].append(
                        merge_base_quality(
                            top['nucleotide'][te],
                            bottom['nucleotide'][be],
                            top['quality'][te],
                            bottom['quality'][be]
                        )
                    )
                t += 1
                b += 1

            merged['nucleotide'] = ''.join(merged['nucleotide'])
            merged['quality'] = ''.join(merged['quality'])
            merged['length'] = len(merged['nucleotide'])
            merged['expected error'] = expected_error(merged['quality'])
            sample['merged'] = merged
            self.report_merge(sample)

    def pick_sample(self, sample):
        def segment_to_fastq(segment):
            return to_fastq(segment['id'], segment['nucleotide'], segment['quality'])

        picked = None
        if self.paired:
            if 'merged' in sample:
                self.trim_segment(sample['merged'])
                picked = sample['merged']

            elif self.node['retain unmerged']:
                self.trim_segment(sample['forward'])
                self.trim_segment(sample['reverse'])
                if  sample['forward']['expected error'] < sample['reverse']['expected error'] and \
                    sample['forward']['length'] >= self.node['minimum read length']:
                    picked = sample['forward']

                elif sample['reverse']['length'] >= self.node['minimum read length']:
                    picked = self.complement_segment(sample['reverse'])
                    # if we use the reverse read we need to directionally trim it again
                    self.trim_segment(picked)
            else:
                self.log.info('dropping unmerged read pair {}\n{}\n{}'.format(sample['id'], segment_to_fastq(sample['forward']), segment_to_fastq(sample['reverse'])))
        else:
            self.trim_segment(sample['forward'])
            picked = deepcopy(sample['forward'])

        # drop reads that are shorter than minimum
        if picked is not None and picked['length'] < self.node['minimum read length']:
            self.log.info('dropping short read {} of length {}\n{}'.format(picked['id'], picked['length'], segment_to_fastq(picked)))
            picked = None

        # drop reads with expected error above the threshold on the minimal length
        if picked is not None:
            effective_error = expected_error(picked['quality'][:self.node['minimum read length']])
            if effective_error > self.node['error threshold']:
                self.log.info('dropping low quality {} read {} with expected error {:.3f}\n{}'.format(picked['origin'], sample['id'], effective_error, segment_to_fastq(picked)))
                picked = None

        if picked is not None:
            picked['id'] = sample['id']
            self.place_in_output(picked)

    def place_in_output(self, segment):
        for term in [
            'effective length',
            'effective start',
            'effective end'
        ]:
            if term in segment:
                del segment[term]
        self.output.append(segment)

    def trim_segment(self, segment):
        length = segment['length']
        width = self.node['trimming window']
        threshold = self.node['trimming error threshold']

        if self.node['trim end']:
            window = segment['quality'][-width:]
            position = 0
            while expected_error(window) > threshold and length >= self.node['minimum read length']:
                length -= 1
                position -= 1
                window = segment['quality'][-width + position:position]

            if position < 0:
                segment['quality'] = segment['quality'][:position]
                segment['nucleotide'] = segment['nucleotide'][:position]

        if self.node['trim start']:
            window = segment['quality'][:width]
            position = 0
            while expected_error(window) > threshold and length >= self.node['minimum read length']:
                length -= 1
                position += 1
                window = segment['quality'][position:width + position]

            if position > 0:
                segment['quality'] = segment['quality'][position:]
                segment['nucleotide'] = segment['nucleotide'][position:]

        segment['length'] = len(segment['nucleotide'])
        segment['expected error'] = expected_error(segment['quality'])

    def score_base_alignment(self, ni, nj, qi, qj):
        score = self.configuration.enumeration['nucleotide substitution'][ni][nj]

        # if qualities were supplied correct the score
        if qi is not None and qj is not None:
            correction =  (phred_to_quality(qi) + phred_to_quality(qj)) / self.node['substitution quality correction']
            score = math.floor(correction * score)

        return score

    def align(self, top, bottom):
        best = None
        for offset in range(top['effective length']):
            score = 0
            mismatch = 0
            t = top['effective start'] + offset
            b = bottom['effective start']
            length = min(top['effective end'] - t, bottom['effective end'] - b)
            if length < self.node['merge overlap threshold']:
                break

            while t < top['effective end'] and b < bottom['effective end']:
                score += self.score_base_alignment(top['nucleotide'][t], bottom['nucleotide'][b], top['quality'][t], bottom['quality'][b])
                if top['nucleotide'][t] != bottom['nucleotide'][b]:
                    mismatch += 1
                t += 1
                b += 1

            if best is None:
                best = {
                    'mismatch': mismatch,
                    'score': score,
                    'offset': offset,
                    'length': length,
                }
            elif score > best['score']:
                best['mismatch'] = mismatch
                best['score'] = score
                best['offset'] = offset
                best['length'] = length

        if best is not None:
            best['top'] = top
            best['bottom'] = bottom

        return best

    def ambiguate(self, i, j):
        return self.ambiguity[i][j]

    def complement_segment(self, segment):
        reversed = {
            'id': segment['id'],
            'nucleotide': ''.join([ self.configuration.enumeration['reverse complement'][n] for n in segment['nucleotide'] ][::-1]),
            'quality': segment['quality'][::-1],
            'length' : segment['length'],
            'expected error': segment['expected error'],
            'origin': segment['origin'],
        }

        if 'effective start' in segment:
            reversed['effective start'] = segment['length'] - segment['effective end']
        if 'effective end' in segment:
            reversed['effective end'] = segment['length'] - segment['effective start']
        if 'effective length' in segment:
            reversed['effective length'] = segment['effective length']
        if 'adapter' in segment:
            reversed['adapter'] = {
                'distance': segment['adapter']['distance'],
                'start': segment['length'] - segment['adapter']['end'],
                'end': segment['length'] - segment['adapter']['start'],
                'length':  segment['adapter']['length']
            }
        return reversed

class SequenceDenoiser(object):
    def __init__(self, pipeline, request=None):
        self.log = logging.getLogger('Noise')
        self.node = {
            'output path': None,
            'iteration': 5,
            'abundance ratio': 10,
            'read length': 133,
            'similarity threshold': 2,
            'ambiguous code': 'BdhKMNRSVWY',
            'statistic': {
                'count': 0
            }
        }
        self.output = None
        self.cache = {}
        if request is not None:
            for term in [
                'limit',
            ]:
                if term in request.instruction:
                    self.node[term] = request.instruction[term]

    def open(self):
        if self.node['output path'] is not None:
            self.output = io.open(self.node['output path'], 'w')
        else:
            self.output = sys.stdout

    def close(self):
        if self.node['output path'] is not None:
            self.output.close()

    @property
    def statistic(self):
        return self.node['statistic']

    @property
    def length(self):
        return self.node['read length']

    def sort(self):
        self.log.info('sorting %s reads in queue', len(self.queue))
        self.queue.sort(key=lambda c: c['error'], reverse=False)
        self.queue.sort(key=lambda c: c['abundance'], reverse=True)

    def add(self, segment):
        if len(segment['nucleotide']) >= self.length:
            self.statistic['count'] += 1
            key = segment['nucleotide'][:self.length]
            abundance = 1
            id = segment['id']

            if ' ' in segment['id']:
                id, comment = segment['id'].split(' ')
                if ':' in comment:
                    abundance, ambiguous = comment.split(':')
                    if abundance is not None:
                        try:
                            abundance = int(abundance)
                        except ValueError:
                            pass

            if key in self.cache:
                cluster = self.cache[key]
                cluster['abundance'] += abundance
            else:
                cluster = {
                    'key': key,
                    'id': id,
                    'abundance': abundance,
                    'ambiguous': False,
                    'quality': [ ],
                    'nucleotide': set(),
                }
                self.cache[key] = cluster
                for c in self.node['ambiguous code']:
                    if c in key:
                        cluster['ambiguous'] = True
                        break
            cluster['quality'].append(segment['quality'])
            cluster['nucleotide'].add(segment['nucleotide'])

    def write(self):
        self.sort()
        for cluster in self.queue:
            id = '{} {}:{}'.format(cluster['id'], cluster['abundance'], 'Y' if cluster['ambiguous'] else 'N')
            self.output.write(to_fastq(id, cluster['key'],cluster['quality']))
            self.output.write('\n')

    def denoise(self):
        self.prepare()
        iteration = 1
        before = len(self.queue) + 1
        while iteration <= self.node['iteration'] and before > len(self.queue):
            before = len(self.queue)
            self.log.info('absorption iteration %d started with %s sequences in queue', iteration, len(self.queue))
            self.bubble()
            iteration += 1
        self.log.info(to_json(self.queue))

    def prepare(self):
        self.queue = list(self.cache.values())
        for cluster in self.queue:
            # consolidate quality
            if cluster['abundance'] > 1:
                quality = []
                count = len(cluster['quality'])
                for i in range(self.length):
                    q = sum([phred_to_quality(j[i]) for j in cluster['quality']])
                    quality.append(quality_to_phred(q / count))
                cluster['quality'] = ''.join(quality)
            else:
                cluster['quality'] = cluster['quality'][0][:self.length]
            cluster['error'] = expected_error(cluster['quality'])

            if cluster['abundance'] == 1:
                cluster['state'] = 'singleton'
            elif cluster['abundance'] >= self.node['abundance ratio']:
                cluster['state'] = 'putative'
            else: 
                cluster['state'] = 'pending'
        self.sort()

    def bubble(self):
        pending = []
        while self.queue and self.queue[-1]['state'] != 'putative':
            cluster = self.queue.pop()
            self.log.debug('pop %s', cluster['id'])
            if self.classify(cluster):
                self.sort()
            else:
                pending.append(cluster)
        self.queue.extend(pending)
        self.sort()

    def report_absorb(self, cluster, putative):
        report = []
        report.append('absorbing {} into {}'.format(cluster['id'], putative['id']))
        report.append('{:<4} {}'.format(putative['abundance'], putative['key']))
        report.append('{:<4} {}'.format('', putative['quality']))
        report.append('{:<4} {}'.format(cluster['abundance'], cluster['key']))
        report.append('{:<4} {}'.format('', cluster['quality']))
        report.append('{:<4} {}'.format('', sequence_diff(putative['key'], cluster['key'])))
        self.log.info('\n'.join(report))

    def absorb(self, cluster, putative):
        self.report_absorb(cluster, putative)

        quality = []
        pa = putative['abundance']
        ca = cluster['abundance']
        for p, c in zip(putative['quality'], cluster['quality']):
            quality.append(quality_to_phred((phred_to_quality(p) * pa + phred_to_quality(c) * ca) / (pa + ca)))
        putative['quality'] = ''.join(quality)

        putative['abundance'] += 1
        putative['nucleotide'] |= cluster['nucleotide']
        if putative['abundance'] >= self.node['abundance ratio']:
            putative['state'] = 'putative'

    def classify(self, cluster):
        result = False
        factor = cluster['abundance'] * self.node['abundance ratio']
        for threshold in range(self.node['similarity threshold']):
            for putative in self.queue:
                if putative['abundance'] < factor:
                    break
                else:
                    distance = hamming(cluster['key'], putative['key'])
                    if putative['abundance'] >= factor and distance <= threshold + 1:
                        self.absorb(cluster, putative)
                        result = True
                        break
            if result:
                break
        return result

class Reference(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Reference')
        self.pipeline = pipeline
        self._lookup = None
        self.node = {
            'key': None,
            'sha1': None,
            'compression': None,
            'remote': None,
            'local': None,
            'path': None,
            'record': None,
        }
        self.node = merge(self.node, node)

    def __str__(self):
        return '{} {}'.format(self.key, self.sha1)

    @property
    def document(self):
        return to_document(self.node)

    @property
    def key(self):
        return self.node['key']

    @property
    def remote(self):
        return self.node['remote']

    @property
    def local(self):
        if not os.path.exists(self.node['local']):
            self.fetch()
        return self.node['local']

    @property
    def path(self):
        return self.node['path']

    @property
    def sha1(self):
        return self.node['sha1']

    @property
    def record(self):
        def append_line_to_record(record, line):
            if record is not None:
                record['body'].write(line.upper())

        def append_record(record):
            if record is not None:
                record['key'] = record['head'].split()[0]
                record['body'].seek(0)
                record['body'] = record['body'].getvalue()
                record['sha1'] = hashlib.sha1(record['body'].encode('utf8')).hexdigest()
                record['sequence'] = Sequence(self.pipeline, { 'nucleotide': record['body'], 'strand': True })
                del record['body']
                self.node['record'][record['key']] = record

        if self.node['record'] is None:
            self.restore()
            if self.node['record'] is None:
                self.node['record'] = {}
                content = self.fetch()
                content = StringIO(self.fetch().decode('utf8'))
                record = None
                for line in content:
                    line = line.strip()
                    if line:
                        if line[0] != '>':
                            append_line_to_record(record, line)
                        else:
                            append_record(record)
                            record = { 
                                'head': line[1:],
                                'body': StringIO(),
                            }
                append_record(record)
                self.persist()
                self.load()
        return self.node['record']

    @property
    def lookup(self):
        if self._lookup is None:
            if self.record is not None:
                self._lookup = { 'key': {}, 'sha1': {} }
                for record in self.record.values():
                    self._lookup['key'][record['key']] = record
                    self._lookup['sha1'][record['sha1']] = record
        return self._lookup

    def record_by_key(self, key):
        if self.lookup is not None and key in self.lookup['key']:
            return self.lookup['key'][key]
        else:
            return None

    def record_by_sha1(self, sha1):
        if self.lookup is not None and sha1 in self.lookup['sha1']:
            return self.lookup['sha1'][sha1]
        else:
            return None

    def load(self):
        if self.node['record'] is not None:
            for record in self.node['record'].values():
                if isinstance(record['sequence'], dict):
                    record['sequence'] = Sequence(self.pipeline, record['sequence'])

    def delete(self):
        if os.path.exists(self.path):
            os.remove(self.path)

    def persist(self):
        try:
            with io.open(self.path, 'wb') as file:
                pickle.dump(self.document, file)
            self.log.debug('persisted reference %s', str(self))
        except OSError as e:
            self.log.error('failed to persist %s to %s', str(self), self.path)
            self.log.debug(str(e))

    def restore(self):
        if os.path.exists(self.path):
            with io.open(self.path, 'rb') as file:
                self.node = pickle.load(file)
                self.load()

    def fetch(self):
        import urllib.request
        import urllib.error
        import http.client

        content = None
        if os.path.exists(self.node['local']):
            with io.open(self.node['local'], 'rb') as file:
                content = file.read()
        else:
            self.log.info('fetching %s', self.remote)
            request = urllib.request.Request(self.remote, None)
            try:
                response = urllib.request.urlopen(request)
            except http.client.BadStatusLine as e:
                self.log.warning('Bad http status error when requesting %s', url)
            except urllib.error.HTTPError as e:
                self.log.warning('Server returned an error when requesting %s: %s', url, e.code)
            except urllib.error.URLError as e:
                self.log.warning('Could not reach server when requesting %s: %s', url, e.reason)
            else:
                content = response.read()
                sha1 = hashlib.sha1(content).hexdigest()
                if sha1 == self.sha1:
                    if len(content) > 3:
                        if content[:3] == b'\x1f\x8b\x08':
                            import gzip
                            self.log.info('decompressing gz archive %s', str(self))
                            content = gzip.decompress(content)

                        if content[:3] == b'\x42\x5a\x68':
                            import bz2
                            self.log.info('decompressing bz2 archive %s', str(self))
                            content = bz2.decompress(content)

                    self.log.info('writing local %s', str(self))
                    with io.open(self.node['local'], 'wb') as file:
                        file.write(content)
                else:
                    raise ValueError('failed checksum verification {}'.format(str(self)))
        return content

class Resolver(object):
    def __init__(self, pipeline):
        self.log = logging.getLogger('Resolver')
        self.pipeline = pipeline
        self._connection = None
        self.cache = {
            'gene':{},
            'sample': {},
            'accession': {},
            'library': {},
            'reference': {},
        }

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def connection(self):
        if self._connection is None:
            try:
                if 'username' in self.configuration.database:
                    self._connection = MongoClient(
                        'mongodb://{}:{}@{}/{}'.format(
                            self.configuration.database['username'],
                            self.configuration.database['password'],
                            self.configuration.database['host'],
                            self.configuration.database['database'],
                        )
                    )
                else:
                    self._connection = MongoClient(
                        'mongodb://{}/{}'.format(
                            self.configuration.database['host'],
                            self.configuration.database['database'],
                        )
                    )
            except pymongo.errors.ConnectionFailure as e:
                self.log.error('failed to establish connection %s', e)
            else:
                self.log.debug('connection established')
        return self._connection

    @property
    def database(self):
        if self.connection is not None:
            return self.connection[self.configuration.database['database']]
        else:
            return None

    def close(self):
        if self._connection is not None:
            self._connection.close()

    def rebuild(self, table):
        objective = []
        existing_collections = self.database.collection_names()
        if table:
            objective = [ t for t in table if t in self.configuration.table ]
        else:
            objective = self.configuration.table.keys()

        for t in objective:
            record = self.configuration.table[t]
            collection = self.database[record['collection']]
            if record['collection'] in existing_collections:
                existing_indexes = collection.index_information()
                for definition in record['index']:
                    if definition['name'] in existing_indexes:
                        self.log.info('dropping index %s on collection %s', definition['name'], record['collection'])
                        collection.drop_index(definition['name'])
                        
            for definition in record['index']:
                self.log.info('building index %s on collection %s', definition['name'], record['collection'])
                collection.create_index([tuple(k) for k in definition['key']], name=definition['name'], unique=definition['unique'])

    def store(self, node):
        if node is not None:
            if 'kind' in node and 'group' in node:
                if 'id' in node:
                    document = {
                        'head': {},
                        'body': node
                    }
                    if node['kind'] == 'gene':
                        for term in [
                            'id',
                            'group',
                            'gene',
                            'strain',
                            'organism',
                            'accession',
                            'framed',
                            'region',
                            'verified',
                            'confirmed',
                            'functionality',
                            'identified',
                        ]:
                            if term in node:
                                document['head'][term] = node[term]

                        gene = Gene(self.pipeline, document)
                        gene.validate()
                        gene.validate_artifact()
                        self.gene_save(gene)

                    elif node['kind'] == 'library':
                        for term in [
                            'id',
                            'strain', 
                            'tissue',
                            'exposure', 
                            'biological repetition', 
                            'technical repetition', 
                        ]:
                            if term in node:
                                document['head'][term] = node[term]

                        library = Library(self.pipeline, document)
                        self.library_save(library)

    def bulk_store_sample(self, buffer):
        bulk = pymongo.bulk.BulkOperationBuilder(self.database['sample'])
        for sample in buffer:
            document = sample.document
            if '_id' in document:
                bulk.find({'_id': document['_id']}).upsert().replace_one(document)
            else:
                bulk.insert(document)
        try:
            bulk.execute()
        except BulkWriteError as e:
            self.log.critical(e.details)
            raise SystemExit()

    def reference_fetch(self, key):
        reference = None
        if key not in self.cache['reference']:
            if key in self.configuration.reference:
                self.cache['reference'][key] = Reference(self, self.configuration.reference[key])

        if key in self.cache['reference']:
            reference = self.cache['reference'][key]

        return reference

    def study_resolve(self, request):
        study = None
        hash = request.hash_for('study')
        name = request.name if request.name else hash

        # explicit name resolution
        if request.name:
            study = self.study_find(request.name)

        # request hash resolution
        if study is None:
            study = self.study_find(hash)

        # if the drop flag was specified and a record exists, drop the record and reset the pointer 
        if study is not None and request.instruction['drop']:
            self.study_drop(study.id)
            study = None

        # if still no study record exists, populate a new one
        if study is None:
            self.log.debug('generating new study %s\n%s', name, to_json(request.document_for('study')))
            study = Study(self.pipeline, None, request)
            cursor = request.cursor_for('sample')
            for node in cursor:
                sample = Sample(self.pipeline, node)
                study.add_sample(sample)
            cursor.close()
            study.finalize()
            self.study_save(study)

        return study

    def study_find(self, keyword):
        study = None
        if keyword is not None:
            document = self.database['study'].find_one({'head.name': keyword})
            if document is None and re.match('^[0-9a-f]{40}$', keyword):
                document = self.database['study'].find_one({'head.id': keyword})
            if document:
                study = Study(self.pipeline, document)
        return study

    def study_drop(self, keyword):
        study = self.study_find(keyword)
        if study:
            study.delete()
            result = self.database['study'].delete_one({ 'head.id': study.id })
            if result:
                self.log.info('dropped study record for %s', study.name)

    def study_save(self, study):
        if study is not None:
            document = to_document({ 'head': study.head })
            existing = self.database['study'].find_one({'head.id': study.id})
            if existing:
                self.log.debug('existing study found for %s', study.name)
                document['_id'] = existing['_id']
            self.database['study'].save(document)
            study.persist()

    def library_save(self, library):
        if library is not None:
            existing = self.library_fetch(library.id)
            if existing:
                self.log.debug('existing library found for %s', library.id)
                library.node['_id'] = existing.node['_id']
                
            self.database['library'].save(library.document)
            if library.id in self.cache['library']:
                del self.cache['library'][library.id]

    def library_fetch(self, id):
        if id not in self.cache['library']:
            document = self.database['library'].find_one({'head.id': id})
            if document:
                library = Library(self.pipeline, document)
                self.cache['library'][id] = library
        if id in self.cache['library']:
            return self.cache['library'][id]
        else:
            return None

    def sample_drop(self, library):
        if library:
            collection = self.database['sample']
            try:
                result = collection.delete_many({ 'head.library': library })
                if result:
                    self.log.info('dropped %d samples from %s', result.deleted_count, library)
            except BulkWriteError as e:
                self.log.critical(e.details)

    def gene_fetch(self, key):
        gene = None
        if key not in self.cache['gene']:
            document = self.database['gene'].find_one({'head.id': key})
            if document:
                self.cache['gene'][key] = Gene(self.pipeline, document)

        if key in self.cache['gene']:
            gene = self.cache['gene'][key]

        return gene

    def gene_save(self, gene):
        if gene is not None:
            existing = self.gene_fetch(gene.id)
            if existing:
                self.log.debug('existing gene found for %s', gene.id)
                gene.node['_id'] = existing.node['_id']
                
            gene.validate()
            gene.validate_artifact()
            self.database['gene'].save(gene.document)
            if gene.id in self.cache['gene']:
                del self.cache['gene'][gene.id]

    def accession_save(self, accession):
        if accession is not None:
            if accession.valid:
                existing = self.database['accession'].find_one({'head.id': accession.id})
                if existing:
                    accession.node['_id'] = existing['_id']
                    self.log.debug('existing accession found for %s', accession.id)
                    
                self.log.debug('saving %s', accession.id)
                self.database['accession'].save(accession.document)
                
                if accession.id in self.cache['accession']:
                    del self.cache['accession'][accession.id]
            else:
                self.log.error('refusing to save invalid accession %s', str(accession))

    def accession_fetch(self, id):
        result = None
        if id not in self.cache['accession']:
            document = self.database['accession'].find_one({'head.id': id})
            if document:
                self.cache['accession'][id] = Accession(self.pipeline, document)
            else:
                ncbi = self.accession_fetch_ncbi(id)
                if ncbi:
                    document = {
                        'head': {
                            'id': ncbi['INSDSeq_primary-accession'],
                            'version': ncbi['INSDSeq_accession-version'],
                            'description': ncbi['INSDSeq_definition'],
                            'simple description': simplify(ncbi['INSDSeq_definition']),
                            'organism name': ncbi['INSDSeq_organism'],
                        },
                        'body': {
                            'sequence': {
                                'nucleotide': ncbi['INSDSeq_sequence'].upper(),
                                'read frame': 0,
                                'strand': True,
                            },
                            'molecule type': ncbi['INSDSeq_moltype']
                        }
                    }
                    
                    if 'INSDSeq_comment' in ncbi:
                        document['body']['comment'] = ncbi['INSDSeq_comment']
                        document['body']['simple comment'] = simplify(ncbi['INSDSeq_comment'])
                        
                    found = False
                    if 'INSDSeq_feature-table' in ncbi:
                        if 'INSDFeature' in ncbi['INSDSeq_feature-table']:
                            for INSDFeature in ncbi['INSDSeq_feature-table']['INSDFeature']:
                                if 'INSDFeature_quals' in INSDFeature:
                                    for INSDFeature_qual in INSDFeature['INSDFeature_quals']:
                                        if 'INSDQualifier' in INSDFeature_qual:
                                            for INSDQualifier in INSDFeature_qual['INSDQualifier']:
                                                if INSDQualifier['INSDQualifier_name'] == 'organism':
                                                    document['head']['organism name'] = INSDQualifier['INSDQualifier_value']
                                                elif INSDQualifier['INSDQualifier_name'] == 'strain':
                                                    document['head']['strain'] = INSDQualifier['INSDQualifier_value']
                                                    found = True
                                                    break
                                        if found: break
                                if found: break
                                
                    if 'strain' not in document['head'] and 'simple description' in document['head']:
                        if 'c57bl/6' in document['head']['simple description']:
                            document['head']['strain'] = 'c57bl6'
                            self.log.info('c57bl6 strain deduced for accession %s from a mention in the description', id)
                            
                    if 'strain' not in document['head'] and 'simple comment' in document['head']:
                        if 'c57bl/6' in document['head']['simple comment']:
                            document['head']['strain'] = 'c57bl6'
                            self.log.info('c57bl6 strain deduced for accession %s from a mention in the description', id)
                            
                    accession = Accession(self.pipeline, document)
                    self.accession_save(accession)
                    self.cache['accession'][id] = accession
            
        if id in self.cache['accession']:
            result = self.cache['accession'][id]
            
        return result

    def accession_fetch_ncbi(self, id):
        import urllib.request
        import urllib.error
        import urllib.parse
        import http.client
        import xmltodict

        def normalize_document(document, transform):
            if isinstance(document, dict):
                transformed = {}
                for k,v in document.items():
                    if k in transform and isinstance(v, dict):
                        transformed[k] = normalize_document([v], transform)
                    else:
                        transformed[k] = normalize_document(v, transform)
                return transformed
                
            elif isinstance(document, list):
                return [ normalize_document(e, transform) for e in document ]
            else:
                return document

        def parse(content):
            transform = [
                "INSDSet",
                "INSDSeq",
                "INSDFeature",
                "INSDFeature_intervals",
                "INSDFeature_quals",
                "INSDQualifier",
                "INSDSeq_references"
            ]
            document = None
            if content:
                node = xmltodict.parse(content.getvalue())
                node = normalize_document(node, transform)
                if 'INSDSet' in node:
                    for INSDSet in node['INSDSet']:
                        if 'INSDSeq' in INSDSet:
                            for INSDSeq in INSDSet['INSDSeq']:
                                if  'INSDSeq_moltype' in INSDSeq and \
                                    INSDSeq['INSDSeq_moltype'] in ['DNA', 'mRNA', 'RNA']:
                                    document = INSDSeq
                                    break
            return document

        def fetch(url):
            content = None
            request = urllib.request.Request(url, None, { 'Accept': 'application/xml' })
            
            try:
                response = urllib.request.urlopen(request)
            except http.client.BadStatusLine as e:
                self.log.warning('Bad http status error when requesting %s', url)
            except urllib.error.HTTPError as e:
                self.log.warning('Server returned an error when requesting %s: %s', url, e.code)
            except urllib.error.URLError as e:
                self.log.warning('Could not reach server when requesting %s: %s', url, e.reason)
            else:
                content = StringIO(response.read().decode('utf8'))
                if content.read(22) == 'Nothing has been found':
                    content = None
                else:
                    content.seek(0)
            return parse(content)

        document = None
        url = self.configuration.expression['ncbi accession url'].format(id)
        node = fetch(url)
        if node and 'INSDSeq_primary-accession' in node:
            document = node
        return document

class Request(object):
    def __init__(self, pipeline, instruction):
        self.log = logging.getLogger('Request')
        self.pipeline = pipeline
        self.node = {
            'instruction': {
                'name': None,
                'kind': None,
                'limit': None,
                'skip': None,
                'action': None,
                'preset': None,
                'profile': None,
            },
            'repertoire': None,
            'hash': None
        }
        if instruction is not None:
            for key,value in instruction.items():
                if value is not None or value not in self.node['instruction']:
                    self.node['instruction'][key] = value

    def reset(self):
        self.node['repertoire'] = None
        self.node['hash'] = None

    @property
    def configuration(self):
        return self.pipeline.configuration

    @property
    def resolver(self):
        return self.pipeline.resolver

    @property
    def instruction(self):
        return self.node['instruction']

    @property
    def action(self):
        return self.instruction['action']

    @property
    def implementation(self):
        implementation = None
        if self.action in self.configuration.interface['implementation']:
            implementation = self.configuration.interface['implementation'][self.action]
        return implementation

    @property
    def format(self):
        if 'format' in self.instruction:
            return self.instruction['format']
        else:
            return None

    @property
    def name(self):
        if self.instruction['name'] is not None:
            return self.instruction['name']
        else:
            return None

    @property
    def kind(self):
        return self.instruction['kind']

    @kind.setter
    def kind(self, value):
        if value != self.instruction['kind']:
            self.instruction['kind'] = value
            self.reset()

    @property
    def limit(self):
        return self.instruction['limit']

    @property
    def skip(self):
        return self.instruction['skip']

    @property
    def repertoire(self):
        if self.node['repertoire'] is None:
            preset = self.preset_for('study')
            if preset is not None:
                self.node['repertoire'] = {}
                template = {
                    'region': 'pivot',
                    'family': 'family',
                    'gene': 'name',
                    'functionality': 'functionality',
                    'reference end': 'end',
                    'reference start': 'start',
                    'reference strand': 'strand',
                    'length': 'length',
                }

                cursor = self.resolver.database['gene'].find (
                    {
                        'head.organism name': preset['organism name'],
                        'head.strain': preset['strain'],
                    }
                )

                for node in cursor:
                    record = { 'key': node['head']['gene'] }
                    for k,v in template.items():
                        if k in node['body']:
                            record[v] = node['body'][k]

                    if record['pivot'] not in self.node['repertoire']:
                        self.node['repertoire'][record['pivot']] = []
                    self.node['repertoire'][record['pivot']].append(record)
        return self.node['repertoire']

    @property
    def profile(self):
        return self.profile_for(self.kind)

    @property
    def query(self):
        return self.query_for(self.kind)

    @property
    def preset(self):
        return self.preset_for(self.kind)

    @property
    def document(self):
        return self.document_for(self.kind)

    @property
    def mongodb_query(self):
        return self.mongodb_query_for(self.kind)

    @property
    def cursor(self):
        return self.cursor_for(self.kind)

    @property
    def count(self):
        return self.count_for(self.kind)

    def profile_for(self, kind):
        profile = None
        if kind is not None:
            if self.instruction['profile'] in self.configuration.profile:
                if kind in self.configuration.profile[self.instruction['profile']]:
                    if self.configuration.profile[self.instruction['profile']][kind] is not None:
                        profile = deepcopy(self.configuration.profile[self.instruction['profile']][kind])
                    else:
                        profile = {}

                    if 'select' not in profile:
                        profile['select'] = {}

                    # potentially override some values from the command line instruction
                    if kind in self.configuration.extension['profile']:
                        if self.configuration.extension['profile'][kind] is not None:
                            if 'select' in self.configuration.extension['profile'][kind]:
                                for key, term in self.configuration.extension['profile'][kind]['select'].items():
                                    if key in self.instruction and self.instruction[key] is not None:
                                        if term['type'] == 'str':
                                            profile['select'][key] = self.instruction[key]

                                        elif term['type'] == 'bool':
                                            if self.instruction[key] == 'Y':
                                                profile['select'][key] = True

                                            elif self.instruction[key] == 'N':
                                                profile['select'][key] = False
        return profile

    def query_for(self, kind):
        query = None
        profile = self.profile_for(kind)
        if profile is not None and 'select' in profile:
            query = {}
            for key,value in profile['select'].items():
                query[key] = value

            # expand library query parameter
            if 'library' in query:
                library = self.resolver.library_fetch(query['library'])
                if library:
                    if library.reference:
                        query[ '__or__' ] = [ ]
                        for reference in library.reference:
                            query[ '__or__' ].append({ 'library': reference })
                        del query['library']
                else:
                    raise ValueError('library {} does not exist'.format(query['library']))
        return query

    def preset_for(self, kind):
        preset = None
        if kind is not None:
            if self.instruction['preset'] in self.configuration.preset:
                if kind in self.configuration.preset[self.instruction['preset']]:
                    preset = deepcopy(self.configuration.preset[self.instruction['preset']][kind])

                    # potentially override some values from the command line instruction
                    if kind in self.configuration.extension['preset']:
                        if self.configuration.extension['preset'][kind] is not None:
                            for key, term in self.configuration.extension['preset'][kind].items():
                                if key in self.instruction and self.instruction[key] is not None:
                                    if term['type'] == 'str':
                                        preset[key] = self.instruction[key]
        return preset

    def document_for(self, kind):
        document = None
        if kind == 'study':
            document = {
                'name': self.name,
                'limit': self.limit,
                'skip': self.skip,
                'query': self.query_for('sample'),
                'preset': self.preset_for('study'),
                'repertoire': self.repertoire,
            }
        else:
            document = {
                'name': self.name,
                'limit': self.limit,
                'skip': self.skip,
                'query': self.query,
                'preset': self.preset,
            }
        return document

    def hash_for(self, kind):
        hash = None
        if kind == 'study':
            o = self.document_for(kind)
            if o['preset'] and 'sort' in o['preset']:
                del o['preset']['sort']  
            o = json.loads(json.dumps(o, sort_keys=True, ensure_ascii=False))
            hash = hashlib.sha1(to_json(o).encode('utf8')).hexdigest()
        return hash

    def mongodb_query_for(self, kind):
        def assemble(node):
            if node is not None:
                if isinstance(node, dict):
                    for key in list(node.keys()):
                        assemble(node[key])
                        if re.match('^__([a-zA-Z]+)__$',key):
                            node['${}'.format(key[2:-2])] = node[key]
                        else:
                            node['head.{}'.format(key)] = node[key]
                        del node[key]
                elif isinstance(node, list):
                    for element in node:
                        assemble(element)
            return node
        return assemble(self.query_for(kind))

    def cursor_for(self, kind):
        cursor = self.resolver.database[kind].find(self.mongodb_query_for(kind))

        if self.limit is not None:
            cursor.limit(self.limit)

        if self.skip is not None:
            cursor.skip(self.skip)

        return cursor

    def count_for(self, kind):
        # is kind and collection the same
        return self.resolver.database[kind].count(self.mongodb_query_for(kind))

class Pipeline(object):
    def __init__(self, configuration):
        self.log = logging.getLogger('Pipeline')
        self.configuration = configuration
        self.resolver = Resolver(self)
        
        for k,v in self.configuration.rss.items():
            heptamer = [ Sequence(self, {'nucleotide': s, 'strand': True, 'read frame': 0}) for s in v['heptamer'] ]
            nonamer = [ Sequence(self, {'nucleotide': s, 'strand': True, 'read frame': 0}) for s in v['nonamer'] ]
            v['3 heptamer'] = re.compile('|'.join([s.nucleotide for s in heptamer]))
            v['3 nonamer'] = re.compile('|'.join([s.nucleotide for s in nonamer]))
            v['5 heptamer'] = re.compile('|'.join([s.reversed.nucleotide for s in heptamer]))
            v['5 nonamer'] = re.compile('|'.join([s.reversed.nucleotide for s in nonamer]))

    def close(self):
        self.resolver.close()

    def execute(self, instruction):
        request = Request(self, instruction)
        if request.implementation:
            action = getattr(self, request.implementation, None)
            if action is not None:
                action(request)
            else:
                raise ValueError('action {} is not implemented'.format(request.action))
        else:
            raise ValueError('action {} is not implemented'.format(request.action))

    def study_collect(self, request):
        if request.format == 'json':
            for key in request.instruction['keys']:
                study = self.resolver.study_find(key)
                study.json(request)

        if request.format == 'csv':
            collection = {
                'study': [],
                'head': [],
                'row': [],
            }

            for key in request.instruction['keys']:
                study = self.resolver.study_find(key)
                if study:
                    collection['study'].append(study)
                    collection['head'].append(study.encode_csv_head(deepcopy(request.preset)))

            if all([collection['head'][0] == h for h in collection['head']]):
                for study in collection['study']:
                    study.encode_csv(request.preset, collection['row'])

                print(to_csv(collection['head'][0]))
                for row in collection['row']:
                    print(to_csv(row))

    def study_expression(self, request):
        pattern = {
            'row': '{index},{pivot},{study},{start},{intensity},{adjusted},{element}',
            'head': {
                'index': 'index',
                'pivot': 'pivot',
                'study': 'study',
                'start': 'start',
                'intensity': 'intensity',
                'adjusted': 'adjusted',
                'element': 'element',
            }
        }
        order = {
            'study': { },
            'pivot': {
                'vh': 3,
                'dh': 2,
                'jh': 1,
            },
        }
        buffer = []

        position = 0
        for name in request.instruction['keys']:
            position += 1
            order['study'][name] = position
            study = self.resolver.study_find(name)
            if study:
                study.tiled_expression(buffer)

        buffer.sort(key=lambda x: -order['pivot'][x['pivot']])
        buffer.sort(key=lambda x: -order['study'][x['study']])
        buffer.sort(key=lambda x: -x['start'])

        index = 0
        current = None
        for row in buffer:
            if current is None or current == row['study']:
                index += 1
                current = row['study']
            row['index'] = index

        print(pattern['row'].format(**pattern['head']))
        for row in buffer:
            print(pattern['row'].format(**row))

    # populate
    def populate(self, request):
        buffer = []
        if request.instruction['path']:
            for path in request.instruction['path']:
                with io.open(path, 'rb') as file:
                    buffer.append(json.loads(file.read().decode('utf8')))
        else:
            buffer.append(json.load(sys.stdin))

        if buffer:
            for document in buffer:
                for kind in [
                    'gene',
                    'library',
                    'study',
                ]:
                    if kind in document:
                        for group in document[kind].keys():
                            for key,record in document[kind][group].items():
                                record['id'] = key
                                record['kind'] = kind
                                record['group'] = group
                                self.resolver.store(record)

    def analyze(self, request):
        count = 0
        if not request.instruction['library']:
            raise ValueError('must specify a library to populate')

        if request.instruction['drop']:
            self.resolver.sample_drop(request.instruction['library'])

        block = Block(self, request)
        while block.fill():
            count += block.size
            # for sample in block.buffer:
            #     if sample.gapped:
            #         print(to_json_document(sample))
            # self.resolver.bulk_store_sample(block.buffer)
            self.log.info('%s so far', count)

    def reanalyze(self, request):
        def flush(buffer, collection):
            if buffer:
                bulk = pymongo.bulk.BulkOperationBuilder(collection)
                for sample in buffer:
                    bulk.find({'_id': sample.document['_id']}).upsert().replace_one(sample.document)
                try:
                    bulk.execute()
                except BulkWriteError as e:
                    self.log.critical(e.details)
                    raise SystemExit()
                return len(buffer)
            else:
                return 0

        cursor = request.cursor_for('sample')
        collection = self.resolver.database['sample']

        count = 0
        buffer = []
        for node in cursor:
            sample = Sample(self, node)
            sample.analyze()
            buffer.append(sample)
            if len(buffer) >= self.configuration.constant['buffer size']:
                count += flush(buffer, collection)
                buffer = []
                self.log.info('%s so far', count)
        cursor.close()
        flush(buffer, collection)

    # drop
    def drop(self, request):
        if request.kind == 'sample':
            self.sample_drop(request)

        if request.kind == 'gene':
            pass

    def sample_drop(self, request):
        self.resolver.sample_drop(request.instruction['library'])

    # count
    def count(self, request):
        if request.kind == 'sample':
            self.sample_count(request)

        if request.kind == 'gene':
            self.gene_count(request)

    def sample_count(self, request):
        print(request.count_for('sample'))

    def gene_count(self, request):
        print(request.count_for('gene'))

    # sample
    def sample(self, request):
        if request.format == 'json':
            self.sample_json(request)

        if request.format == 'fasta':
            self.sample_fasta(request)

        if request.format == 'fastq':
            self.sample_fastq(request)

        if request.format == 'diagram':
            self.sample_alignment(request)

    def sample_json(self, request):
        cursor = request.cursor_for('sample')
        for node in cursor:
            sample = Sample(self, node)
            sample.info()
        cursor.close()

    def sample_fasta(self, request):
        cursor = request.cursor_for('sample')
        for node in cursor:
            sample = Sample(self, node)
            print(sample.fasta)
        cursor.close()

    def sample_fastq(self, request):
        cursor = request.cursor_for('sample')
        for node in cursor:
            sample = Sample(self, node)
            print(sample.fastq)
        cursor.close()

    def sample_alignment(self, request):
        cursor = request.cursor_for('sample')
        for node in cursor:
            sample = Sample(self, node)
            sample.view(request)
        cursor.close()

    def sample_filter(self, request):
        reader = FastqFilteringReader(self, request)
        reader.open()

        while reader.fill() and not reader.enough:
            reader.write()

        reader.close()
        self.log.debug(to_json(reader.statistic))

    def sample_denoise(self, request):
        request.instruction['filtering'] = False
        reader = FastqFilteringReader(self, request)
        denoiser = SequenceDenoiser(self, request)

        denoiser.open()
        reader.open()

        while reader.fill() and not reader.enough:
            for sample in reader.output:
                denoiser.add(sample)

        reader.close()
        self.log.debug(to_json(reader.statistic))
        reader = None

        denoiser.denoise()
        denoiser.write()
        denoiser.close()

    # study
    def study(self, request):
        if request.format == 'csv':
            self.study_csv(request)

        if request.format == 'json':
            self.study_json(request)

        if request.format == 'markdown':
            self.study_markdown(request)

    def study_csv(self, request):
        study = self.resolver.study_resolve(request)
        study.csv(request)

    def study_json(self, request):
        study = self.resolver.study_resolve(request)
        study.json(request)

    def study_markdown(self, request):
        cursor = request.cursor_for('study')
        for node in cursor:
            study = Study(self, node)
            print(str(study))
        cursor.close()

    # library
    def library(self, request):
        if request.format == 'json':
            self.library_json(request)

        if request.format == 'markdown':
            self.library_markdown(request)

    def library_json(self, request):
        cursor = request.cursor_for('library')
        buffer = []
        for node in cursor:
            buffer.append(node)
        cursor.close()
        buffer.sort(key=lambda x: '' if 'technical' not in x else x['technical'])
        buffer.sort(key=lambda x: '' if 'biological' not in x else x['biological'])
        buffer.sort(key=lambda x: '' if 'tissue' not in x['head'] else x['head']['tissue'])
        buffer.sort(key=lambda x: '' if 'exposure' not in x['head'] else x['head']['exposure'])
        print(to_json(buffer))

    def library_markdown(self, request):
        buffer = {}
        order = []
        cursor = request.cursor_for('library')
        cursor.sort (
            [
                ('head.exposure', pymongo.DESCENDING),
                ('head.tissue', pymongo.DESCENDING),
                ('body.biological repetition', pymongo.ASCENDING),
                ('body.technical repetition', pymongo.ASCENDING),
            ]
        )

        for node in cursor:
            library = Library(self, node)
            if library.reference is None:
                if library.count is None or request.instruction['drop']:
                    r = Request(self, {'library': library.id, 'profile': 'all', 'kind': 'sample'})
                    library.count = r.count

                    r = Request(self, {'library': library.id, 'profile': 'p', 'kind': 'sample'})
                    library.productive_count = r.count

                    r = Request(self, {'library': library.id, 'profile': 'valid', 'kind': 'sample'})
                    library.valid_count = r.count

                    self.resolver.library_save(library)

            buffer[library.id] = library
            order.append(library)
        cursor.close()

        for library in buffer.values():
            if library.reference is not None:
                if library.count is None or request.instruction['drop']:
                    library.count = 0
                    library.productive_count = 0
                    library.valid_count = 0
                    for reference in library.reference:
                        if reference in buffer:
                            library.count += buffer[reference].count
                            library.productive_count += buffer[reference].productive_count
                            library.valid_count += buffer[reference].valid_count
                    self.resolver.library_save(library)


        title = ' | '.join (
            ['{:<9}', '{:<3}', '{:<2}', '{:<2}', '{:<30}', '{:12}', '{:8}', '{:8}', '{:<25}']
        ).format (
            'Strain', 'Exp', 'B', 'T', 'Tissue', 'Count', 'Valid', 'Prod', 'ID'
        )
        print('| ' + title)

        separator = '-|-'.join (
            ['{:-<9}', '{:-<3}', '{:-<2}', '{:-<2}', '{:-<30}', '{:-<12}', '{:-<8}', '{:-<8}', '{:-<25}']
        ).format (
            '', '', '', '', '', '', '', '', ''
        )
        print('| ' + separator)
        for library in order:
            print(str(library))

    # gene
    def index(self, request):
        self.index_bwa(request)
        self.index_igblast(request)

    def index_bwa(self, request):
        chain = self.configuration.selected['chain']
        verify_directory(chain['bwa directory'])
        for region in self.configuration.selected['region'].values():
            if region['type'] == 'gene segment':
                instruction = deepcopy(request.instruction)
                instruction['action'] = 'gene'
                instruction['format'] = 'fasta'
                instruction['kind'] = 'gene'
                instruction['region'] = region['key']

                r = Request(self, instruction)
                cursor = r.cursor
                buffer = []
                for node in cursor:
                    buffer.append(node)
                cursor.close()
                buffer.sort(key=lambda x: '' if 'allele' not in x['head'] else x['head']['allele'])
                buffer.sort(key=lambda x: '' if 'gene' not in x['head'] else x['head']['gene'])
                buffer.sort(key=lambda x: '' if 'family' not in x['head'] else x['head']['family'])

                self.log.debug('writing %s region bwa fasta file', region['key'])
                with io.open(region['bwa fasta path'], 'wb') as file:
                    for node in buffer:
                        gene = Gene(self, node)
                        file.write(gene.to_fasta(request.instruction['flanking'], self.configuration.constant['fasta line length']).encode('utf8'))
                        file.write(b'\n')

                self.log.debug('building %s region bwt index', region['key'])
                command = Command(self, 'bwa index')
                command.args.append(region['bwa fasta path'])
                output, error = command.process.communicate()

    def index_igblast(self, request):
        chain = self.configuration.selected['chain']
        verify_directory(chain['igblast database directory'])
        # copy the internal data directory into the igblast database directory
        command = Command(self, 'rsync')
        command.parameter['--recursive'] = None
        command.args.append(self.configuration.constant['igblastn internal data'])
        command.args.append(chain['igblast directory'])
        command.process.communicate()

        # build blast db files
        for region in self.configuration.selected['region'].values():
            if region['type'] == 'gene segment':
                instruction = deepcopy(request.instruction)
                instruction['action'] = 'gene'
                instruction['format'] = 'fasta'
                instruction['kind'] = 'gene'
                instruction['region'] = region['key']

                r = Request(self, instruction)
                cursor = r.cursor
                buffer = []
                for node in cursor:
                    buffer.append(node)
                cursor.close()
                buffer.sort(key=lambda x: '' if 'allele' not in x['head'] else x['head']['allele'])
                buffer.sort(key=lambda x: '' if 'gene' not in x['head'] else x['head']['gene'])
                buffer.sort(key=lambda x: '' if 'family' not in x['head'] else x['head']['family'])

                self.log.debug('writing %s region igblast fasta file', region['key'])
                with io.open(region['igblast fasta path'], 'wb') as file:
                    for node in buffer:
                        gene = Gene(self, node)
                        file.write(gene.to_fasta(request.instruction['flanking'], self.configuration.constant['fasta line length']).encode('utf8'))
                        file.write(b'\n')

                self.log.debug('building %s region blast index', region['key'])
                command = Command(self, 'makeblastdb', region['makeblastdb'])
                command.cwd = chain['igblast database directory']
                command.process.communicate()

    def gene(self, request):
        if request.format == 'json':
            self.gene_json(request)

        if request.format == 'html':
            self.gene_html(request)

        if request.format == 'fasta':
            self.gene_fasta(request)

        if request.format == 'markdown':
            self.gene_markdown(request)

    def gene_json(self, request):
        cursor = request.cursor_for('gene')
        buffer = []
        for node in cursor:
            document = node['body'].copy()
            document.update(node['head'])
            buffer.append(document)
        cursor.close()
        buffer.sort(key=lambda x: '' if 'allele' not in x else x['allele'])
        buffer.sort(key=lambda x: '' if 'gene' not in x else x['gene'])
        buffer.sort(key=lambda x: '' if 'family' not in x else x['family'])
        buffer.sort(key=lambda x: '' if 'strain' not in x else x['strain'])
        print(to_json(buffer))

    def gene_html(self, request):
        cursor = request.cursor_for('gene')
        buffer = []
        for node in cursor:
            buffer.append(Gene(self, node))
        cursor.close()
        buffer.sort(key=lambda x: '' if 'allele' not in x.body else x.body['allele'])
        buffer.sort(key=lambda x: '' if 'gene' not in x.body else x.body['gene'])
        buffer.sort(key=lambda x: '' if 'family' not in x.body else x.body['family'])
        buffer.sort(key=lambda x: '' if 'strain' not in x.body else x.body['strain'])
        buffer.sort(key=lambda x: '' if 'reference start' not in x.body else x.body['reference start'])

        print("""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" >
            <html xmlns="http://www.w3.org/1999/xhtml">
            <head>
            <meta content="text/html; charset=UTF-8" http-equiv="Content-Type" />""")

        if request.instruction['title'] is not None: print('<title>{}</title>'.format(request.instruction['title']))

        print("""<style type="text/css">
                html {
                    font-family: Courier;
                    font-size: 12px;
                    background-color: #ffffff;
                    color: #000000;
                }
                body {
                    margin: 1em;
                }
                * {
                    padding: 0;
                    margin: 0;
                }
                img, object, a {
                    border-width: 0;
                    outline: none;
                }
                .wrap {
                    margin: 0 auto;
                }
                .section {
                    white-space: nowrap;
                    padding: 0.5em;
                    margin: 0.5em 0;
                }
                .genename {
                    font-weight: bold;
                    color: #222;
                }
                .gene {
                    display: inline;
                    color: #555;
                }
                .codon {
                    color: #5A7251;
                }
                .codon .cdr1,
                .codon .cdr2,
                .codon .cdr3 {
                    color: #3A5231;
                }
                .codon .fr1,
                .codon .fr2,
                .codon .fr3 {
                    color: #7A9271;
                }
                .frame {
                    font-weight: bold;
                }
                .nonamer {
                    color: #FBAC21;
                }
                .spacer {
                    color: #9ECA2B;
                }
                .heptamer {
                    color: #82A5D6;
                }
                .gap {
                    color: #F14329;
                }
                .nonamer,
                .spacer,
                .heptamer,
                .gap,
                .preamble,
                .fr1,
                .fr2,
                .fr3,
                .cdr1,
                .cdr2,
                .cdr3 {
                    display: inline;
                }
                .gene .cdr1,
                .gene .cdr2,
                .gene .cdr3 {
                    color: #333;
                }
                .gene .fr1,
                .gene .fr2,
                .gene .fr3 {
                    color: #AAA;
                }
                .implied {}
                .BN000872 {}
                .placeholder {
                    text-decoration: line-through;
                }
                .pattern {
                    text-decoration: underline;
                }
                a {
                    text-decoration: none;
                }
                a:link,a:visited {
                    color: #2d4855;
                }
                a:hover {
                    color: #0f70c5;
                }
            </style>
            </head>
            <body><div class="wrap">""")

        for gene in buffer: gene.html()
        print("""</div></body></html>""")

    def gene_fasta(self, request):
        cursor = request.cursor_for('gene')
        buffer = []
        for node in cursor:
            buffer.append(node)
        cursor.close()
        buffer.sort(key=lambda x: '' if 'allele' not in x['body'] else x['body']['allele'])
        buffer.sort(key=lambda x: '' if 'gene' not in x['body'] else x['body']['gene'])
        buffer.sort(key=lambda x: '' if 'family' not in x['body'] else x['body']['family'])
        buffer.sort(key=lambda x: '' if 'strain' not in x['body'] else x['body']['strain'])
        for node in buffer:
            gene = Gene(self, node)
            print(gene.to_fasta(request.instruction['flanking'], self.configuration.constant['fasta line length']))
        cursor.close()

    def gene_rss(self, request):
        cursor = request.cursor_for('gene')
        for node in cursor:
            gene = Gene(self, node)
            gene.check_rss(request.instruction['flanking'], request.instruction['distance'])
            self.resolver.gene_save(gene)
        cursor.close()

    def gene_align(self, request):
        instruction = {
            'record': {}, 
            'total': 0, 
            'flank': request.instruction['flanking'],
        }
        
        cursor = request.cursor_for('gene')
        for node in cursor:
            gene = Gene(self, node)
            instruction['record'][gene.id] = { 'gene': gene }
            instruction['total'] += 1
        cursor.close()
        
        # execute blat and add the hits to each gene sequence in the buffer
        self.log.info('aligning %s sequences with %s flanking', instruction['total'], instruction['flank'])
        blat = Blat(self, instruction)
        blat.search()
        for k,query in blat.instruction['record'].items():
            gene = query['gene']
            if 'summary' in query and query['summary']:
                if len(query['summary']) == 1:
                    hit = query['summary'][0]
                    gene.assign_reference_alignment(hit)
                    self.resolver.gene_save(gene)
                else:
                    self.log.warning('gene %s aligned to multiple locations', gene.id)
            else:
                self.log.error('no satisfactory alignment found for gene %s',  gene.id)

    def gene_markdown(self, request):
        cursor = request.cursor_for('gene')
        cursor.sort ([('body.reference start', pymongo.DESCENDING)])
        buffer = []
        for node in cursor:
            gene = Gene(self, node)
            o = {
                'gene': gene.body['gene'],
                'functionality': gene.body['functionality'],
                'length': gene.body['length'],
                'family': gene.body['family'],
                'reference start': gene.body['reference start'],
                'gap': None,
            }
            buffer.append(o)
        cursor.close()

        for index, gene in enumerate(buffer):
            if index > 0:
                gene['gap'] = buffer[index - 1]['reference start'] - gene['reference start']

        title = ' | '.join ( ['{:<13}', '{:<12}', '{:<1}', '{:<3}', '{:<12}', '{:<6}'] ).format ( 'gene', 'family', 'F', 'len', 'ref', 'gap' )
        print('| ' + title)

        separator = '-|-'.join (['{:-<13}', '{:-<12}', '{:-<1}', '{:-<3}', '{:-<12}', '{:-<6}']).format ('', '', '', '', '', '')
        print('| ' + separator)

        for gene in buffer:
            print('| ' + ' | '.join(['{:<13}', '{:<12}', '{:<1}', '{:<3}', '{:<12}', '{:<6}']).format(
                str(gene['gene']), 
                str(gene['family']), 
                str(gene['functionality']), 
                str(gene['length']), 
                str(gene['reference start']), 
                '' if gene['gap'] is None else str(gene['gap']), 
            ))

    # accession
    def accession(self, request):
        if request.format == 'json':
            self.accession_json(request)

        if request.format == 'fasta':
            self.accession_fasta(request)

    def accession_fasta(self, request):
        cursor = request.cursor_for('accession')
        for node in cursor:
            accession = Accession(self, node)
            print(accession.fasta)
        cursor.close()

    def accession_json(self, request):
        cursor = request.cursor_for('accession')
        for node in cursor:
            accession = Accession(self, node)
            print(accession.to_json())
        cursor.close()

    # db
    def rebuild(self, request):
        self.resolver.rebuild(request.instruction['table'])

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    
    configuration = Configuration()
    command = CommandLineParser(configuration.interface)
    if command.sectioned and command.action is None:
        command.help()

    else:
        configuration.select(command.instruction)
        logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])
        pipeline = Pipeline(configuration)
        try:
            pipeline.execute(command.instruction)

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
