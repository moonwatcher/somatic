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

import xmltodict
import urllib.request
import urllib.parse
import urllib.error
import urllib.parse
import urllib.request
import urllib.error
import http.client


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

def to_simple_json(node):
    def handler(o):
        result = None
        if isinstance(o, datetime):
            result = o.isoformat()
        if isinstance(o, ObjectId):
            result = str(o)
        return result
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4, default=handler)

def parse_match(match):
    if match is not None:
        return dict(((k.replace('_', ' '),v) for k,v in match.groupdict().items() if k and v))
    else:
        return None

def to_fastq(id, sequence, quality):
    if id:
        if sequence:
            if quality:
                return '{}\n{}\n+\n{}'.format(id, sequence, quality)
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

def transform_to_document(node):
    if isinstance(node, set):
        return [ transform_to_document(v) for v in node ]
    if isinstance(node, list):
        return [ transform_to_document(v) for v in node ]
    elif isinstance(node, dict):
        return dict([ (k,transform_to_document(v)) for k,v in node.items() ])
    elif isinstance(node, Sample):
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
    elif isinstance(node, numpy.ndarray):
        if node.dtype == numpy.float:
            return [ float(v) for v in node ]
        elif node.dtype == numpy.int:
            return [ int(v) for v in node ]
        else:
            return node
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
    document = transform_to_document(node)
    return json.dumps(document, sort_keys=True, ensure_ascii=False, indent=4, default=handler)

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

def load_configuration():
    def merge(node, other):
        result = other
        if node is not None:
            result = node
            if other is not None:
                if isinstance(node, dict):
                    if isinstance(other, dict):
                        for k,v in other.items():
                            if k in result:
                                result[k] = merge(result[k], v)
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

    def check(element):
        result = None
        if isinstance(element, dict):
            if 'enabled' not in element or element['enabled']:
                if 'enabled' in element:
                    del element['enabled']

                for k in list(element.keys()):
                    element[k] = check(element[k])
                result = element

        elif isinstance(element, list):
            result = []
            for o in element:
                checked = check(o)
                if checked is not None:
                    result.append(checked)
        else:
            result = element

        return result

    setting = None
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'configuration/setting.json')
    with io.open(path, 'rb') as file:
        setting = json.loads(file.read().decode('utf8'))

    configuration = None
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'configuration/core.json')
    with io.open(path, 'rb') as file:
        configuration = json.loads(file.read().decode('utf8'))
        configuration['cache'] = setting['cache']
        configuration['database'] = setting['database']
        configuration['reference'] = setting['reference']
        for k,v in setting['constant'].items():
           configuration['constant'][k] = v

    configuration['command']['blat']['arguments'][0] = configuration['constant']['blat binary']
    configuration['command']['blat']['arguments'][1] = configuration['reference']['mouse chromosome 12']

    configuration['command']['igblast']['arguments'][0] = configuration['constant']['igblastn binary']
    configuration['command']['igblast']['cwd'] = configuration['constant']['igblastn database']

    if configuration is not None:
        configuration['interface']['prototype']['profile']['parameter']['choices'] = list(configuration['profile'].keys())
        configuration['interface']['prototype']['profile']['parameter']['help'] = 'one of: {}'.format(', '.join(sorted(configuration['profile'].keys())))

        if 'study' in configuration:
            for key in list(configuration['study']['feature'].keys()):
                feature = deepcopy(configuration['study']['default']['feature'])
                for column in list(feature['column'].keys()):
                    feature['column'][column] = merge(deepcopy(configuration['study']['default']['column']), feature['column'][column])
                    feature['column'][column]['key'] = column

                feature['key'] = key
                feature = merge(feature, configuration['study']['feature'][key])
                configuration['study']['feature'][key] = feature

            for name in list(configuration['study']['preset'].keys()):
                preset = deepcopy(configuration['study']['default']['preset'])
                preset['key'] = name
                preset = merge(preset, configuration['study']['preset'][name])

                for key in list(preset['feature'].keys()):
                    feature = deepcopy(configuration['study']['feature'][key])
                    feature = merge(feature, preset['feature'][key])
                    preset['feature'][key] = feature

                row = []
                for feature_key in preset['order']:
                    if feature_key in preset['feature']:
                        feature = { 'key': feature_key, 'column': [] }
                        for column_key in preset['feature'][feature_key]['order']:
                            if column_key in preset['feature'][feature_key]['column']:
                                feature['column'].append(preset['feature'][feature_key]['column'][column_key])
                        row.append(feature)
                preset['row'] = check(row)
                del preset['feature']
                del preset['order']
                configuration['study']['preset'][name] = preset

        configuration['interface']['prototype']['preset']['parameter']['choices'] = list(configuration['study']['preset'].keys())
        configuration['interface']['prototype']['preset']['parameter']['help'] = 'one of {}'.format(', '.join(sorted(configuration['study']['preset'].keys())))

        configuration['expression'] = {
            'ncbi accession url': 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&dopt=gbc_xml&val={}',
            'blat hit': re.compile(
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
            ),
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
            'nucleotide sequence': re.compile('^[ACGTRYKMSWBDHVN]+$', re.IGNORECASE)
        }
        configuration['diagram'] = {
            'prototype': {
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
        }

    # load a reverse complement table for ambiguity code
    configuration['complement'] = {}
    for k,v in configuration['iupac nucleic acid notation'].items():
        configuration['complement'][k] = v['reverse']

    # load collection of codons possibly encoding a stop codon
    configuration['stop codon repertoire'] = possibly_stop_codon(configuration)

    configuration['interface']['implementation'] = {}
    for action in configuration['interface']['section']['action']:
        configuration['interface']['implementation'][action['instruction']['name']] = action['implementation']

    return configuration

def possibly_stop_codon(configuration):
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
    
    seed = [ k for k,v in configuration['nucleic to amino'].items() if v == '*' ]
    feature = {
        'codon': dict([ (c, { 'motif': [] }) for c in seed ]),
        'space': set()
    }

    for triplet, codon in feature['codon'].items():
        for nucleotide in triplet:
            motif = []
            for k,n in configuration['iupac nucleic acid notation'].items():
                if nucleotide in n['option']:
                    motif.append(k)
            codon['motif'].append(motif)

        codon['possible'] = expand(codon['motif'])
        feature['space'] |= set(codon['possible'])
    return feature['space']


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


class Blat(object):
    def __init__(self, pipeline, instruction):
        self.log = logging.getLogger('Blat')
        self.pipeline = pipeline
        self.instruction = instruction
        self.target = StringIO()
        self.fasta = []

        for record in self.instruction['record'].values():
            record['flanking'] = record['gene'].flanking_accession_query(instruction['flank'])
            self.fasta.append(
                to_fasta(
                    record['flanking']['id'],
                    record['flanking']['sequence'].nucleotide,
                    None,
                    self.configuration['constant']['fasta line length']
                )
            )
        self.fasta = '\n'.join(self.fasta)

        with io.open(self.instruction['target path'], 'r') as file:
            for line in file:
                line = line.strip()
                if line[0] != '>':
                    self.target.write(line.upper())
        self.target.seek(0)

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

    def crop(self, query):
        if 'hit' in query:
            for hit in query['hit']:
                orientation = self.orientation(query, hit['flanked'])
                    
                # extract both query and target sequence for the contiguous blocks
                for block in hit['flanked']['block']:
                    self.target.seek(block['target start'])
                    block['target sequence'] = self.target.read(block['size'])
                    block['query sequence'] = orientation['flanking'].nucleotide[block['query start']:block['query start'] + block['size']]
                    
                hit['objective'] = {
                    'block': [],
                    'score': 0,
                    'identical': 0.0,
                    'match': 0,
                    'mismatch': 0,
                    'block count': 0,
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

                        self.target.seek(o['target start'])
                        o['target sequence'] = self.target.read(o['size'])
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
                self.target.seek(position['target'])
                gap = block['target start'] - position['target']
                hit['target sequence'] += self.target.read(gap)
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
        for i in range(flanked['block count']):
            block = {
                'size': flanked['block size'][i],
                'query start': flanked['query block start'][i],
                'target start': flanked['target block start'][i],
            }
            self.target.seek(block['target start'])
            block['target sequence'] = self.target.read(block['size'])
            block['query sequence'] = orientation['flanking'].nucleotide[block['query start']:block['query start'] + block['size']]
            flanked['block'].append(block)

        del flanked['block size']
        del flanked['query block start']
        del flanked['target block start']
        del flanked['query name']

        if 'hit' not in query: query['hit'] = []
        query['hit'].append({'flanked': self.normalize_hit(flanked, orientation)})

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

    def search(self):
        command = self.configuration['command']['blat']
        process = Popen(
            args=command['arguments'],
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate(input=self.fasta.encode('utf8'))
        if output:
            buffer = StringIO(output.decode('utf8'))
            for line in buffer:
                hit = parse_match(self.configuration['expression']['blat hit'].search(line.strip()))
                if hit:
                    query = self.instruction['record'][hit['query name']]
                    self.parse_flanked_hit(query, hit)

            # for every query sort the results by the objective score and than the flanked score
            for query in self.instruction['record'].values():
                if 'hit' in query:
                    self.crop(query)
                    query['hit'] = sorted(query['hit'], key=lambda x: x['flanked']['identical'], reverse=True)
                    query['hit'] = sorted(query['hit'], key=lambda x: x['objective']['identical'], reverse=True)
                    query['hit'] = sorted(query['hit'], key=lambda x: x['flanked']['score'], reverse=True)
                    query['hit'] = sorted(query['hit'], key=lambda x: x['objective']['score'], reverse=True)
                    self.summarize(query)


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
        document = transform_to_document(self.node)
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
        if start is not None:
            start = max(start, 0)
            end = self.length if end is None else min(end, self.length)
            if end > start:
                cropped = {
                    'nucleotide': self.nucleotide[start:end],
                    'read frame': (3 - (start - self.read_frame)%3)%3,
                    'strand': self.strand
                }
                if 'quality' in self.node:
                    cropped['quality'] = self.quality[start:end]
                sequence = Sequence(self.pipeline, cropped)
        return sequence

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
            return ''

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
            return ''

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
    def classified(self):
        classified = None
        if self.codon:
            classified = []
            for codon in self.codon:
                acid = self.configuration['iupac amino acid notation'][codon]
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
        return len(self.nucleotide)

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
        return transform_to_document(self.node)

    @property
    def fasta(self):
        return self.to_fasta(self.configuration['constant']['fasta line length'])

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
        return transform_to_document(self.node)

    @property
    def json(self):
        return to_json(self.node)

    @property
    def fasta(self):
        return to_fasta(self.id, self.sequence.nucleotide, None, self.configuration['constant']['fasta line length'])

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
            reference = self.pipeline.resolver.reference_fetch('mouse chromosome 12')
            if reference is not None:
                query = {}
                query['id'] = self.id
                query['reference start'] = self.body['reference start']
                query['reference end'] = self.body['reference end']
                query['reference strand'] = self.body['reference strand']
                query['flank start'] = max(query['reference start'] - flank, 0)
                query['flank end'] = min(query['reference end'] + flank, reference.length)
                query['flank length'] = query['flank end'] - query['flank start']
                if self.body['reference strand']:
                    query['start'] = query['reference start'] - query['flank start']
                    query['end'] = query['reference end'] - query['flank start']
                    query['sequence'] = reference.crop(query['flank start'], query['flank end'])
                else:
                    query['start'] = query['flank end'] - query['reference end']
                    query['end'] = query['flank end'] - query['reference start']
                    query['sequence'] = reference.crop(query['flank start'], query['flank end']).reversed
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
                    if artifact != 'preamble' and fragment.read_frame != 0:
                        self.log.info('artifact %s for %s is out of frame %s', artifact, self.id, fragment.read_frame)
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
            if self.region == 'VH':
                self._check_rss_3(flanking, 'VH', distance)

            elif self.region == 'DH':
                self._check_rss_5(flanking, 'DH', distance)
                self._check_rss_3(flanking, 'DH', distance)

            elif self.region == 'JH':
                self._check_rss_5(flanking, 'JH', distance)

    def html(self):
        buffer = []
        buffer.append('<div class="section">')
        buffer.append('<div class="genename">{} {}</div>'.format(self.id, self.functionality))

        self._codon_html(buffer)        
        if self.region == 'VH':
            self._nucleotide_html(buffer)
            self._html_3(buffer, 'VH')

        elif self.region == 'DH':
            self._html_5(buffer, 'DH')
            self._nucleotide_html(buffer)
            self._html_3(buffer, 'DH')

        elif self.region == 'JH':
            self._html_5(buffer, 'JH')
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
            if self.region == 'VH' and framed_artifact:
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
        if self.region == 'VH':
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
        pattern = self.configuration['rss'][region][name]
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
                    'id': id,
                    'library': library,
                    'valid': True
                }
            }
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

    def make_default_hit(self):
        return {
            'valid': True,
            'framed': False,
            'in frame': False,
            'gapped': False,
            'picked': False,
            'region': None,
            'subject id': None,
            'subject strand': self.strand,
            'query start': None,
            'query end': None,
            'query': None,
        }

    def add_igblast_hit(self, hit):
        if hit is not None:
            if 'subject id' in hit:
                gene = self.pipeline.resolver.gene_fetch(hit['subject id'])
                # switch the hit to a zero based coordinate system
                if gene:
                    hit['query start'] -= 1
                    if hit['subject strand'] == gene.sequence.strand:
                        hit['subject start'] -= 1
                    else:
                        hit['subject start'] = gene.length - hit['subject start']
                        hit['subject end'] = gene.length - hit['subject end'] + 1
                else:
                    self.invalidate_hit(hit, 'unknown gene')
            self._add_hit(hit)

    def analyze(self, strain=None):
        self._reset()
        self._load_gene()
        self._pick_jh_region(strain)
        self._pick_vh_region(strain)
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
        self.head['CDR3 length'] = None
        self.head['DH length'] = None
        self.head['JH length'] = None
        self.head['VH length'] = None
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

        if 'hit' in self.body:
            self.body['hit'] = [ hit for hit in self.body['hit'] if hit['region'] in self.configuration['region'].keys()]

        for hit in self.hit:
            hit['valid'] = True
            hit['picked'] = False
            hit['framed'] = False
            hit['in frame'] = False
            hit['score'] = hit['bit score']
            if 'gap openings' in hit and hit['gap openings'] > 0:
                hit['gapped'] = True
                self.head['gapped'] = True
                self.invalidate_hit(hit, 'gapped')
            else:
                hit['gapped'] = False
                hit['query'] = self.sequence.crop(hit['query start'], hit['query end'])

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
        if 'uuid' not in hit:
            hit['uuid'] = str(uuid.uuid4())

        if hit['uuid'] not in self.index:
            self.index[hit['uuid']] = hit

        for k in [
            'query',
            'subject',
            '3 chew',
            '5 chew'
        ]:
            if k in hit and isinstance(hit[k], dict):
                hit[k] = Sequence(self.pipeline, hit[k])

    def _load_gene(self):
        for hit in self.hit:
            if hit['valid'] and 'subject id' in hit:
                gene = self.pipeline.resolver.gene_fetch(hit['subject id'])
                if gene:
                    hit['framed'] = gene.framed
                    hit['functionality'] = gene.functionality
                    hit['strain'] = gene.strain
                    hit['gene'] = gene.body['gene']
                    hit['family'] = gene.body['family']
                    hit['allele'] = gene.body['allele']
                    if hit['subject strand'] == gene.sequence.strand:
                        hit['subject'] = gene.sequence.crop(hit['subject start'], hit['subject end'])
                        hit['reference strand'] = gene.reference_strand
                    else:
                        hit['subject'] = gene.sequence.reversed.crop(hit['subject start'], hit['subject end'])
                        hit['reference strand'] = not gene.reference_strand

                    if hit['reference strand']:
                        hit['reference start'] = gene.reference_start + hit['subject start']
                        hit['reference end'] = gene.reference_start + hit['subject end']
                    else:
                        hit['reference start'] = gene.reference_start + gene.length - hit['subject end']
                        hit['reference end'] = gene.reference_start + gene.length - hit['subject start']

                    # only DH regions are allowed to align to the opposite strand
                    if hit['region'] != 'DH' and hit['subject strand'] != gene.sequence.strand:
                        self.invalidate_hit(hit, 'wrong strand')
                else:
                    self.invalidate_hit(hit, 'unknown gene')

    def _frame_hits(self):
        if self.framed:
            for h in self.hit:
                if h['valid']:
                    if not h['framed']:
                        h['query'].read_frame = 2 - (h['query start'] - self.sequence.read_frame - 1) % 3
                        h['framed'] = True

                    if not h['in frame'] and 'subject' in h:
                        if h['query'].read_frame == h['subject'].read_frame:
                            h['in frame'] = True

    def _pick_hit(self, hit):
        def pick_region_primary(hit):
            if hit['region'] not in self.body['primary']:
                self._primary = None
                hit['primary'] = True
                self.body['primary'][hit['region']] = hit['uuid']

        def pick_framing_region(hit):
            # first hit to be picked sets the frame for the sample
            if (hit is not None and \
                not self.framed and \
                hit['framed'] and \
                'query start' in hit and \
                'subject' in hit and \
                'subject strand' in hit):
                self.head['framed'] = True
                self.body['framed by'] = hit['uuid']
                self.sequence.read_frame = (hit['query start'] + hit['subject'].read_frame) % 3
                self.head['strand'] = hit['subject strand']
                # self.log.debug('frame set for %s by %s', self.id, hit['subject id'])

        if hit is not None and hit['query'] is not None:
            hit['picked'] = True
            pick_region_primary(hit)
            pick_framing_region(hit)
            self._frame_hits()

    def _pick_region(self, region, candidate, strain=None):
        picked = False
        if candidate:
            top = None
            search = { 'valid': [] }
            for hit in candidate:
                if (hit['valid'] and
                    hit['region'] == self.configuration['region'][region]['name'] and
                    hit['identical'] >= self.configuration['region'][region]['minimum identity'] and
                    hit['alignment length'] >= self.configuration['region'][region]['minimum alignment']):
                    search['valid'].append(hit)
                    
            if search['valid']:
                top = search['valid']
                top.sort(reverse=True, key=lambda x: x['score'])
                top = [ hit for hit in top if hit['score'] == top[0]['score'] ]
                
                if len(top) > 1:
                    search['framed'] = []
                    for hit in top:
                        frame = self.sequence.read_frame if self.framed else hit['subject'].read_frame
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

                # if there are multiple matches and at least one is functional, keep only the functional
                if len(top) > 1:
                    search['functional'] = [ hit for hit in top if hit['functionality'] == 'F' ]
                    if search['functional']:
                        top = search['functional']
            if top:
                # prefer a functional gene for a primary F > O > P
                top.sort(key=lambda x: x['functionality'])
                for hit in top: self._pick_hit(hit)
                picked = True
        return picked

    def _pick_jh_region(self, strain):
        return self._pick_region('JH', self.hit, strain)

    def _pick_vh_region(self, strain):
        return self._pick_region('VH', self.hit, strain)

    def _pick_dh_region(self):
        jh = None if 'JH' not in self.primary else self.primary['JH']
        vh = None if 'VH' not in self.primary else self.primary['VH']

        if vh is not None and jh is not None:
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
                        gene = self.pipeline.resolver.gene_fetch(hit['subject id'])
                        trimmed = { 'overlap': 0 }
                        for k,v in hit.items():
                            if k not in ('query', 'subject', 'overlap'):
                                trimmed[k] = v

                        # check for overlap with VH
                        if vh['query end'] > hit['query start']:
                            trimmed['query start'] = vh['query end']
                            overlap = trimmed['query start'] - hit['query start']
                            trimmed['subject start'] += overlap
                            trimmed['overlap'] += overlap
                            
                        # check for overlap with JH
                        if hit['query end'] > jh['query start']:
                            trimmed['query end'] = jh['query start']
                            overlap = hit['query end'] - trimmed['query end']
                            trimmed['subject end'] -= overlap
                            trimmed['overlap'] += overlap
                            
                        # correct the score for the overlap and trim the sequences
                        trimmed['score'] = hit['bit score'] - trimmed['overlap'] * self.configuration['region']['DH']['overlap penalty factor']
                        trimmed['alignment length'] = trimmed['query end'] - trimmed['query start']
                        trimmed['query'] = self.sequence.crop(trimmed['query start'], trimmed['query end'])
                        
                        if trimmed['subject strand'] == gene.sequence.strand:
                            trimmed['subject'] = gene.sequence.crop(trimmed['subject start'], trimmed['subject end'])
                        else:
                            trimmed['subject'] = gene.sequence.reversed.crop(trimmed['subject start'], trimmed['subject end'])
                            
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
                            if longest >= self.configuration['region']['DH']['minimum alignment']:
                                search['trimmed'].append(original)
                    if search['trimmed']:
                        top = search['trimmed']
            self._pick_region('DH', top)

    def _check_v_j_framing(self):
        jh = None if 'JH' not in self.primary else self.primary['JH']
        vh = None if 'VH' not in self.primary else self.primary['VH']
        if jh is not None and vh is not None:
            if jh['in frame'] and vh['in frame']:
                self.head['in frame'] = True
        else:
            self.invalidate('could not establish a VH JH pair')

    def _identify_junction(self, left, right, name):
        def identify_palindrome(junction, left, right):
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
                self.head['p count'] += p
                self.head['n count'] += n
                if p > 0: self.head['palindromic'] = True
                if n > 0:
                    if p > 0:
                        junction['palindrome ratio'] = float(n) / float(p)
                    else:
                        junction['palindrome ratio'] = 1.0
                else:
                    junction['palindrome ratio'] = 0.0

        if left is not None and right is not None:
            if left['query end'] < right['query start']:
                junction = self.make_default_hit()
                junction['region'] = name
                junction['subject id'] = name
                junction['query start'] = left['query end']
                junction['query end'] = right['query start']
                junction['query'] = self.sequence.crop(left['query end'], right['query start'])
                identify_palindrome(junction, left['query'], right['query'])
                self._add_hit(junction)
                self._pick_hit(junction)

    def _identify_v_d_junction(self):
        vh = None if 'VH' not in self.primary else self.primary['VH']
        dh = None if 'DH' not in self.primary else self.primary['DH']
        self._identify_junction(vh, dh, 'V-D')

    def _identify_d_j_junction(self):
        dh = None if 'DH' not in self.primary else self.primary['DH']
        jh = None if 'JH' not in self.primary else self.primary['JH']
        self._identify_junction(dh, jh, 'D-J')

    def _identify_v_j_junction(self):
        if 'DH' not in self.primary:
            vh = None if 'VH' not in self.primary else self.primary['VH']
            jh = None if 'JH' not in self.primary else self.primary['JH']
            self._identify_junction(vh, jh, 'V-J')

    def _check_for_stop_codon(self):
        if self.effective.codon and self.framed:
            if '*' in self.effective.codon:
                self.head['premature'] = True

    def _identify_cdr3(self):
        def locate_cdr3_start(vh):
            position = None
            gene = self.pipeline.resolver.gene_fetch(vh['subject id'])
            if 'cdr3' in gene.body:
                gene_cdr3 = gene.body['cdr3']
                if gene_cdr3['reference strand']:
                    offset = gene_cdr3['reference start'] - vh['reference start']
                else:
                    offset = vh['reference end'] - gene_cdr3['reference end']
                position = vh['query start'] + offset
                
                # because the gene annotation considers the CDR3 to start after the Cycteine
                position -= 3
                
                if position < 0 or position > self.sequence.length - 3:
                    self.make_comment('CDR3 start position out of range {}'.format(position))
                    position = None
                else:
                    acid = self.sequence.nucleotide[position:position+3]
                    if acid in self.configuration['nucleic to amino'] and self.configuration['nucleic to amino'][acid] == 'C':
                        self.head['conserved cdr3 c'] = True
                    else:
                        self.make_comment('framing cycteine not conserved {}'.format(acid))
            return position

        def locate_cdr3_end(jh):
            position = None
            gene = self.pipeline.resolver.gene_fetch(jh['subject id'])
            if 'cdr3' in gene.body:
                gene_cdr3 = gene.body['cdr3']
                if gene_cdr3['reference strand']:
                    offset = gene_cdr3['reference end'] - jh['reference start']
                else:
                    offset = jh['reference end'] - gene_cdr3['reference start']
                position = jh['query start'] + offset

                # because the gene annotation considers the CDR3 to end before the Tryptophan
                position += 3

                if position < 0 or position > self.sequence.length - 3:
                    self.make_comment('CDR3 end position out of range {}'.format(position))
                    position = None

                acid = self.sequence.nucleotide[position-3:position]
                if acid in self.configuration['nucleic to amino'] and self.configuration['nucleic to amino'][acid] == 'W':
                    self.head['conserved cdr3 w'] = True
                else:
                    self.make_comment('framing tryptophan not conserved {}'.format(acid))
            return position

        jh = None if 'JH' not in self.primary else self.primary['JH']
        vh = None if 'VH' not in self.primary else self.primary['VH']
        if vh is not None and jh is not None:
            start = locate_cdr3_start(vh)
            end = locate_cdr3_end(jh)
            query = self.sequence.crop(start, end)
            if query:
                cdr3 = self.make_default_hit()
                cdr3['region'] = 'CDR3'
                cdr3['subject id'] = 'CDR3'
                cdr3['query start'] = start
                cdr3['query end'] = end
                cdr3['query'] = query
                self._add_hit(cdr3)
                self._pick_hit(cdr3)

                cdr3['charge'] = 0
                cdr3['weight'] = 0

                if cdr3['framed'] and cdr3['query'].codon:
                    for acid in cdr3['query'].codon:
                        amino = self.configuration['iupac amino acid notation']
                        if 'charge' in amino[acid]:
                            cdr3['charge'] += amino[acid]['charge']
                        if 'weight' in amino[acid]:
                            cdr3['weight'] += amino[acid]['weight']

    def _identify_chewback(self):
        for hit in self.hit:
            if hit['valid'] and hit['region'] in self.configuration['region'].keys():
                gene = self.pipeline.resolver.gene_fetch(hit['subject id'])
                if gene:
                    if hit['subject strand'] == gene.sequence.strand:
                        ref = gene.sequence
                    else:
                        ref = gene.sequence.reversed

                    if hit['region'] == 'VH' or hit['region'] == 'DH':
                        if hit['subject end'] < ref.length:
                            hit['3 chew'] = ref.crop(hit['subject end'], ref.length)
                            
                    if hit['region'] == 'JH' or hit['region'] == 'DH':
                        if hit['subject start'] > 0:
                            hit['5 chew'] = ref.crop(0, hit['subject start'])

    def _check_for_productivity(self):
        jh = None if 'JH' not in self.primary else self.primary['JH']
        vh = None if 'VH' not in self.primary else self.primary['VH']
        cdr3 = None if 'CDR3' not in self.primary else self.primary['CDR3']
        if (vh is not None and \
            jh is not None and \
            cdr3 is not None and \
            self.in_frame and \
            not self.premature and \
            self.conserved_cdr3_cycteine and \
            self.conserved_cdr3_tryptophan):
            if vh['functionality'] == 'F':
                if jh['functionality'] == 'F':
                    self.head['productive'] = True
                else:
                    self.make_comment('leading JH {} is non functional {}'.format(jh['allele'], jh['functionality']))
            else:
                self.make_comment('leading VH {} is non functional {}'.format(vh['allele'], vh['functionality']))

    def _analyze_quality(self):
        for name, region in self.primary.items():
            self.head['{} length'.format(name)] = region['query'].length
            region['average phread'] = float(sum(region['query'].phred)) / float((region['query'].length))

        self.head['average phred'] = float(sum(self.effective.phred)) / float((self.effective.length))

    def invalidate(self, message):
        self.head['valid'] = False
        self.make_comment(message)

    def make_comment(self, message):
        if 'comment' not in self.body:
            self.body['comment'] = []

        self.body['comment'].append(message)
        self.log.info('sample %s : %s', self.id, message)

    def make_hit_comment(self, hit, message):
        if 'comment' not in hit:
            hit['comment'] = []
        hit['comment'].append(message)
        self.log.info('%dbp %s hit to %s on %s : %s', hit['alignment length'], hit['region'], hit['subject id'], self.id, message)

    def invalidate_hit(self, hit, message):
        hit['valid'] = False
        self.make_hit_comment(hit, message)

    def reverse(self):
        self.body['sequence'] = self.body['sequence'].reversed

    def view(self, profile):
        diagram = Diagram(self.pipeline, self, profile)
        print(diagram.draw())

    def info(self, profile):
        print(self.json)

    @property
    def has_dh(self):
        return 'DH' in self.primary

    @property
    def observable(self):
        return 'CDR3' in self.primary and \
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
                        'id': self.trimmed_id,
                        'cdr3 length': float(self.primary['CDR3']['query'].length),
                        'cdr3 charge': float(self.primary['CDR3']['charge']),
                        'cdr3 weight': float(self.primary['CDR3']['weight']),
                    }
                }

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
                    breakdown['region'] = { 'VH': [], 'DH': [], 'JH': [] }
                    for hit in self.hit:
                        if hit['picked'] and hit['region'] in breakdown['region'].keys():
                            breakdown['region'][hit['region']].append(hit)
                    breakdown['combination'] = len(breakdown['region']['VH']) * len(breakdown['region']['DH']) * len(breakdown['region']['JH'])
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
                    chew = [ v['3 chew'].length for v in breakdown['region']['VH'] if '3 chew' in v ]
                    if chew: breakdown['feature']['v3 chew'] = float(numpy.mean(chew))
                    chew = [ j['5 chew'].length for j in breakdown['region']['JH'] if '5 chew' in j ]
                    if chew: breakdown['feature']['j5 chew'] = float(numpy.mean(chew))
                    chew = [ d['3 chew'].length for d in breakdown['region']['DH'] if '3 chew' in d ]
                    if chew: breakdown['feature']['d3 chew'] = float(numpy.mean(chew))
                    chew = [ d['5 chew'].length for d in breakdown['region']['DH'] if '5 chew' in d ]
                    if chew: breakdown['feature']['d5 chew'] = float(numpy.mean(chew))
                    breakdown['feature']['chew'] = \
                        breakdown['feature']['v3 chew'] + \
                        breakdown['feature']['j5 chew'] + \
                        breakdown['feature']['d3 chew'] + \
                        breakdown['feature']['d5 chew']

                    # junction
                    if 'V-D' in self.primary:
                        breakdown['feature']['vd length'] = float(self.primary['V-D']['query'].length)
                        breakdown['feature']['vd n count'] = float(self.primary['V-D']['palindrome'].count('N'))
                        breakdown['feature']['vd p count'] = float(self.primary['V-D']['palindrome'].count('P'))
                        breakdown['feature']['n count'] += breakdown['feature']['vd n count']
                        breakdown['feature']['p count'] += breakdown['feature']['vd p count']
                    if 'D-J' in self.primary:
                        breakdown['feature']['dj length'] = float(self.primary['D-J']['query'].length)
                        breakdown['feature']['dj n count'] = float(self.primary['D-J']['palindrome'].count('N'))
                        breakdown['feature']['dj p count'] = float(self.primary['D-J']['palindrome'].count('P'))
                        breakdown['feature']['n count'] += breakdown['feature']['dj n count']
                        breakdown['feature']['p count'] += breakdown['feature']['dj p count']

                    for vh in breakdown['region']['VH']:
                        for dh in breakdown['region']['DH']:
                            for jh in breakdown['region']['JH']:
                                observation = deepcopy(breakdown['feature'])
                                observation['VH'] = vh['gene']
                                observation['DH'] = dh['gene']
                                observation['JH'] = jh['gene']
                                observation['vh family'] = vh['family']
                                observation['dh family'] = dh['family']
                                self._observation.append(observation)
                else:
                    breakdown['region'] = { 'VH': [], 'JH': [] }
                    for hit in self.hit:
                        if hit['picked'] and hit['region'] in breakdown['region'].keys():
                            breakdown['region'][hit['region']].append(hit)
                    breakdown['combination'] = len(breakdown['region']['VH']) * len(breakdown['region']['JH'])
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
                    chew = [ float(j['5 chew'].length) for j in breakdown['region']['JH'] if '5 chew' in j ]
                    if chew: breakdown['feature']['j5 chew'] = float(numpy.mean(chew))
                    chew = [ float(v['3 chew'].length) for v in breakdown['region']['VH'] if '3 chew' in v ]
                    if chew: breakdown['feature']['v3 chew'] = float(numpy.mean(chew))
                    breakdown['feature']['chew'] = \
                        breakdown['feature']['j5 chew'] + \
                        breakdown['feature']['v3 chew']

                    # junction
                    if 'V-J' in self.primary:
                        breakdown['feature']['vj length'] = float(self.primary['V-J']['query'].length)
                        breakdown['feature']['vj n count'] = float(self.primary['V-J']['palindrome'].count('N'))
                        breakdown['feature']['vj p count'] = float(self.primary['V-J']['palindrome'].count('P'))

                    for vh in breakdown['region']['VH']:
                        for jh in breakdown['region']['JH']:
                            observation = deepcopy(breakdown['feature'])
                            observation['VH'] = vh['gene']
                            observation['JH'] = jh['gene']
                            observation['vh family'] = vh['family']
                            self._observation.append(observation)

        return self._observation

    @property
    def cdr3_sequence(self):
        if 'CDR3' in self.primary:
            return self.primary['CDR3']
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
        return transform_to_document(self.node)

    @property
    def json(self):
        return to_json(self.node)

    @property
    def fasta(self):
        return to_fasta(self.id, self.sequence.nucleotide, None, self.configuration['constant']['fasta line length'])

    @property
    def fastq(self):
        return to_fastq(self.id, self.sequence.nucleotide, self.sequence.quality)

    @property
    def id(self):
        return self.head['id']

    @property
    def trimmed_id(self):
        return self.head['id'].split()[0]

    @property
    def length(self):
        return self.sequence.length

    @property
    def strand(self):
        return self.head['strand']

    @property
    def library(self):
        return self.head['library']

    @property
    def sequence(self):
        return self.body['sequence']

    @property
    def effective(self):
        if 'effective' not in self.body:
            start = self.sequence.length
            end = 0
            for hit in self.hit:
                if hit['valid']:
                    start = min(start, hit['query start'])
                    end = max(end, hit['query end'])

            if start > end:
                self.body['effective'] = self.sequence
                self.body['effective query start'] = 0
                self.body['effective query end'] = self.sequence.length
            else:
                self.body['effective query start'] = start
                self.body['effective query end'] = end
                self.body['effective'] = self.sequence.crop(start, end)
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
    def empty(self):
        return not self.count > 0

    @property
    def document(self):
        document = transform_to_document(self.node)
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

    @property
    def mean(self):
        return self.node['residual']

    def add_observation(self, observation):
        self.count += 1
        self.weight.append(observation['weight'])
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

    def collect(self, template={}, buffer=[]):
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
            template['weight'] = numpy.sum(self.weight)
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
            if 'body' in self.node:
                if 'repertoire' in self.body:
                    for key in list(self.repertoire.keys()):
                        self.body['repertoire'][key] = dict([(x['key'], x) for x in self.repertoire[key]])

                if 'pivot' in self.body and self.body['pivot'] is not None:
                    self.body['pivot'] = Pivot(self, self.body['pivot'], None, self.rotate)

        elif request is not None and request.preset:
            self.node = {
                'head' :{
                    'id': request.hash,
                    'name': request.name,
                },
                'body': {
                    'request': request.document,
                }
            }
            self.reset()
        else:
            raise ValueError('insufficient information to create study')

    def delete(self):
        path = os.path.join(self.configuration['cache']['path'], self.id)
        if os.path.exists(path):
            os.remove(path)

    def persist(self):
        path = os.path.join(self.configuration['cache']['path'], self.id)
        try:
            verify_directory(self.configuration['cache']['path'])
            with io.open(path, 'wb') as file:
                pickle.dump(self.document, file)
            self.log.debug('pickled study %s to %s', self.id, path)
        except OSError as e:
            self.log.error('failed writing pickled study %s to %s', self.id, path)
            self.log.debug(str(e))

    def restore(self):
        path = os.path.join(self.configuration['cache']['path'], self.id)
        if os.path.exists(path):
            with io.open(path, 'rb') as file:
                self.node = pickle.load(file)
                self.load()

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
        return self.body['request']

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
            self.body['pivot'] = Pivot(self, None, None, self.rotate)

        if self.body['pivot'] is None:
            self.restore()

        return self.body['pivot']

    @property
    def repertoire(self):
        return self.body['repertoire']

    @property
    def document(self):
        node = transform_to_document(self.node)
        if 'repertoire' in node['body']:
            for key in list(node['body']['repertoire'].keys()):
                node['body']['repertoire'][key] = list(node['body']['repertoire'][key].values())
        return node

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

    def collect(self, template, buffer=[]):
        buffer = self.pivot.collect({ 'depth': template['depth'] })
        return buffer

    def table(self, preset=None):
        def construct(template):
            row = self.rotate[:template['depth']] + [ 'weight' , 'portion' ]
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

        def sort(template, structure, buffer):
            for sort in template['sort']:
                if sort['column'] in structure['name']:
                    index = structure['name'][sort['column']]
                buffer = sorted(buffer, key=lambda x: x[index], reverse=sort['reverse'])
            return buffer

        def table(template, collection, buffer=[]):
            weight = float(self.count)
            for record in collection:
                row = record['rotate'] + [ record['weight'], record['weight'] / weight ]
                for feature in template['row']:
                    for column in feature['column']:
                        if feature['key'] in record['residual'] and column['key'] in record['residual'][feature['key']]:
                            row.append(record['residual'][feature['key']][column['key']])
                        else:
                            row.append(None)
                buffer.append(row)
            return buffer

        def stringify(template, buffer):
            result = []
            for row in buffer:
                result.append([str(column) for column in row])
            return result

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
        collection = self.collect(template)

        buffer = table(template, collection)
        buffer = sort(template, structure, buffer)
        buffer = stringify(template, buffer)

        print(','.join(structure['head']))
        for row in buffer:
            print(','.join(row))

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
            tile['intensity'] += predicate['weight']
            tile['pivot'].append(predicate)
        return intensity

    def tiled_expression(self, buffer):
        density = {
            'VH': {
                'order': 3,
                'interval': 70000,
            },
            'DH': {
                'order': 2,
                'interval': 2500,
            },
            'JH': {
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
        print(to_simple_json(transform_to_document(self.node)))


class Diagram(object):
    def __init__(self, pipeline, sample, profile):
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
        if profile in self.configuration['profile'] and 'diagram' in self.configuration['profile'][profile]:
            p = self.configuration['profile'][profile]['diagram']
        else:
            p = self.configuration['profile']['default']['diagram']
            
        for k,v in p['track'].items():
            self.query[k] = v
            
        for track in sample.hit:
            if track['region'] not in self.configuration['region'].keys()  or all([k in track and track[k] == v for k,v in self.query.items()]):
                self.track.append(track)
                
        for k in p['feature']:
            if k in self.configuration['diagram']['prototype']:
                feature = dict(self.configuration['diagram']['prototype'][k])
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
        b.append(self.sample.id)
        if 'VH' in self.sample.primary:
            vh = self.sample.primary['VH']
            b.append('{} : {}'.format(vh['region'], vh['gene']))

        if 'DH' in self.sample.primary:
            dh = self.sample.primary['DH']
            b.append('{} : {}'.format(dh['region'], dh['gene']))

        if 'JH' in self.sample.primary:
            jh = self.sample.primary['JH']
            b.append('{} : {}'.format(jh['region'], jh['gene']))

        if 'CDR3' in self.sample.primary:
            cdr3 = self.sample.primary['CDR3']
            if 'charge' in cdr3 and 'weight' in cdr3:
                b.append('{} : {:.3} {}'.format(cdr3['region'], cdr3['charge'], cdr3['weight']))

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
            if offset > 0:
                buffer.write(' ' * offset)
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

    def to_fasta(self, query=None):
        buffer = StringIO()
        for sample in self.buffer:
            if query is None or all([k in sample.head and sample.head[k] == v for k,v in query.items()]):
                buffer.write(to_fasta(sample.id, sample.sequence.nucleotide, None, self.configuration['constant']['fasta line length']))
                buffer.write('\n')
        buffer.seek(0)
        return buffer

    def parse_igblast(self, buffer):
        def parse_igblast_hit(line):
            hit = None
            match = self.configuration['expression']['igblast hit'].search(line)
            if match:
                hit = parse_match(match)
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
                    
            elif state == 2 and line.startswith(self.configuration['expression']['igblast reversed query']):
                # this means the hits will be for the reverse complement strand
                sample.reverse()
                
            elif state == 2 and line.startswith('# Hit table '):
                # this means the hit table is on the next line
                state = 3
                
            elif state == 3 and not line.startswith('#'):
                hit = parse_igblast_hit(line)
                if hit is not None:
                    sample.add_igblast_hit(hit)

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


class Library(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Library')
        self.pipeline = pipeline
        self.node = node

    def __str__(self):
        return '| ' + ' | '.join(['{:<9}', '{:<3}', '{:<2}', '{:<2}', '{:<30}', '{:12}', '{:8}', '{:8}', '{:<25}', '{}']).format(
            '' if 'strain' not in self.head else str(self.head['strain']), 
            '' if 'exposure' not in self.head else str(self.head['exposure']), 
            '' if 'biological repetition' not in self.body else str(self.body['biological repetition']), 
            '' if 'technical repetition' not in self.body else str(self.body['technical repetition']), 
            '' if 'tissue' not in self.head else str(self.head['tissue']),
            str(self.count), 
            '0' if self.count == 0 else '{:.4}'.format(self.valid_count/self.count * 100.0), 
            '0' if self.valid_count == 0 else '{:.4}'.format(self.productive_count / self.valid_count * 100.0), 
            self.id, 
            '' if self.path is None else self.path
        )

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
    def path(self):
        if 'path' in self.body:
            return self.body['path']
        else:
            return None

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

    @property
    def document(self):
        return transform_to_document(self.node)


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
                self._connection = MongoClient(
                    'mongodb://{}:{}@{}/{}'.format(
                        self.configuration['database']['username'],
                        self.configuration['database']['password'],
                        self.configuration['database']['host'],
                        self.configuration['database']['database'],
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
            return self.connection[self.configuration['database']['database']]
        else:
            return None

    def close(self):
        if self._connection is not None:
            self._connection.close()

    def rebuild(self, table):
        objective = []
        existing_collections = self.database.collection_names()
        if table:
            objective = [ t for t in table if t in self.configuration['table'] ]
        else:
            objective = self.configuration['table'].keys()

        for t in objective:
            record = self.configuration['table'][t]
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

    def reference_fetch(self, name):
        reference = None
        if name in self.configuration['reference']:
            if name not in self.cache['reference']:
                record = StringIO()
                with io.open(self.configuration['reference'][name], 'r') as file:
                    for line in file:
                        line = line.strip()
                        if line[0] != '>':
                            record.write(line.upper())
                record.seek(0)
                sequence = {
                    'nucleotide': record.getvalue(),
                    'read frame': 0,
                    'strand': True
                }
                self.cache['reference'][name] = Sequence(self.pipeline, sequence)
            if name in self.cache['reference']:
                reference = self.cache['reference'][name]
        else:
            self.log.error('unknown reference %s', name)
        return reference

    def study_resolve(self, request):
        request.kind = 'sample'
        study = None

        # explicit name resolution
        if request.instruction['name'] is not None:
            study = self.study_find(request.instruction['name'])

        # request hash resolution
        if study is None:
            study = self.study_find(request.hash)

        # if the drop flag was specified and a record exists, drop the record and reset the pointer 
        if study is not None and request.instruction['drop']:
            self.study_drop(study.id)
            study = None

        # if still no study record exists, populate a new one
        if study is None:
            study = Study(self.pipeline, None, request)
            self.log.debug('populating study %s', study.name)
            cursor = request.cursor('analyzed_sample')
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
            document = study.document
            document['body']['pivot'] = None
            existing = self.database['study'].find_one({'head.id': study.id})
            if existing:
                self.log.debug('existing study found for %s', study.name)
                document['_id'] = existing['_id']
            self.database['study'].save(document)
            study.persist()

    def library_store(self, node):
        if node is not None:
            library = Library(self.pipeline, node)
            self.library_save(library)

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

    def analyzed_sample_drop(self, library):
        if library:
            collection = self.database['analyzed_sample']
            try:
                result = collection.delete_many({ 'head.library': library })
                if result:
                    self.log.info('dropped %d samples from %s', result.deleted_count, library)
            except BulkWriteError as e:
                self.log.critical(e.details)

    def make_repertoire(self, request):
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

        organism = request['preset']['organism name']
        strain = request['preset']['strain']
        repertoire = {}
        query = {
            'head.organism name': organism,
            'head.strain': strain
        }
        cursor = self.database['gene'].find(query)
        for node in cursor:
            record = { 'key': node['head']['gene'] }
            for k,v in template.items():
                if k in node['body']:
                    record[v] = node['body'][k]

            if record['pivot'] not in repertoire:
                repertoire[record['pivot']] = []
            repertoire[record['pivot']].append(record)
        return repertoire

    def gene_fetch(self, id):
        if id not in self.cache['gene']:
            document = self.database['gene'].find_one({'head.id': id})
            if document:
                gene = Gene(self.pipeline, document)
                self.cache['gene'][id] = gene
        if id in self.cache['gene']:
            return self.cache['gene'][id]
        else:
            return None

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

    def gene_store(self, node):
        if node is not None:
            document = { 'head': {}, 'body': node }
            for k in [
                'accession',
                'framed',
                'region',
                'id',
                'verified',
                'confirmed',
                'functionality',
                'organism name',
                'strain',
                'identified',
                'gene',
            ]:
                if k in node:
                    document['head'][k] = node[k]

            gene = Gene(self.pipeline, document)
            gene.validate()
            gene.validate_artifact()
            self.gene_save(gene)

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
                            document['head']['strain'] = 'C57BL/6'
                            self.log.info('C57BL/6 strain deduced for accession %s from a mention in the description', id)
                            
                    if 'strain' not in document['head'] and 'simple comment' in document['head']:
                        if 'c57bl/6' in document['head']['simple comment']:
                            document['head']['strain'] = 'C57BL/6'
                            self.log.info('C57BL/6 strain deduced for accession %s from a mention in the description', id)
                            
                    accession = Accession(self.pipeline, document)
                    self.accession_save(accession)
                    self.cache['accession'][id] = accession
            
        if id in self.cache['accession']:
            result = self.cache['accession'][id]
            
        return result

    def accession_fetch_ncbi(self, id):
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
                log.warning('Bad http status error when requesting %s', url)
            except urllib.error.HTTPError as e:
                log.warning('Server returned an error when requesting %s: %s', url, e.code)
            except urllib.error.URLError as e:
                log.warning('Could not reach server when requesting %s: %s', url, e.reason)
            else:
                content = StringIO(response.read().decode('utf8'))
                if content.read(22) == 'Nothing has been found':
                    content = None
                else:
                    content.seek(0)
            return parse(content)

        document = None
        url = self.configuration['expression']['ncbi accession url'].format(id)
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
            'query': None,
            'profile': None,
            'preset': None,
            'repertoire': None,
            'hash': None
        }
        if instruction is not None:
            for key,value in instruction.items():
                if value is not None or value not in self.node['instruction']:
                    self.node['instruction'][key] = value

    def reset(self):
        self.node['profile'] = None
        self.node['query'] = None
        self.node['preset'] = None
        self.node['repertoire'] = None
        self.node['hash'] = None

    @property
    def document(self):
        return {
            'name': self.name,
            'skip': self.skip,
            'limit': self.limit,
            'query': self.query,
            'preset': self.preset,
            'repertoire': self.repertoire,
        }

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
    def query(self):
        if self.node['query'] is None:
            self.node['query'] = {}

            # load a query profile
            if self.profile is not None:
                for key,value in self.profile.items():
                    self.node['query'][key] = value

            # numeric and string instruction parameters that go into the query
            for key in [
                'id',
                'names',
                'region',
                'format',
                'strain',
                'library',
                'functionality',
            ]:
                if key in self.instruction and self.instruction[key] is not None:
                    self.node['query'][key] = self.instruction[key]

            # boolean instruction parameters that go into the query
            for key in [
                'valid',
                'gapped',
                'framed',
                'in frame',
                'premature',
                'productive',
            ]:
                if key in self.instruction and self.instruction[key] is not None:
                    if self.instruction[key] == 'Y':
                        self.node['query'][key] = True

                    elif self.instruction[key] == 'N':
                        self.node['query'][key] = False

            # expand library query parameter
            if 'library' in self.node['query']:
                library = self.resolver.library_fetch(self.node['query']['library'])
                if library:
                    if library.reference:
                        self.node['query'][ '__or__' ] = [ ]
                        for reference in library.reference:
                            self.node['query'][ '__or__' ].append({ 'library': reference })
                        del self.node['query']['library']
                else:
                    raise ValueError('library {} does not exist'.format(self.node['query']['library']))

        return self.node['query']

    @property
    def hash(self):
        if self.node['hash'] is None:
            o = {
                'limit': self.limit,
                'skip': self.skip,
                'query': self.query,
                'preset': deepcopy(self.preset),
                'repertoire': self.repertoire,
            }

            if o['preset'] and 'sort' in o['preset']:
                del o['preset']['sort']  

            o = json.loads(json.dumps(o, sort_keys=True, ensure_ascii=False))
            self.node['hash'] = hashlib.sha1(to_json(o).encode('utf8')).hexdigest()
        return self.node['hash']

    @property
    def name(self):
        if self.instruction['name'] is not None:
            return self.instruction['name']
        else:
            return self.hash

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
    def profile(self):
        if self.node['profile'] is None and self.kind is not None:
            if self.instruction['profile'] in self.configuration['profile']:
                if self.kind in self.configuration['profile'][self.instruction['profile']]:
                    self.node['profile'] = deepcopy(self.configuration['profile'][self.instruction['profile']][self.kind])
                else:
                    raise ValueError('profile {} is undefined for {}'.format(self.instruction['profile'], self.kind))
            else:
                raise ValueError('profile {} is undefined'.format(self.instruction['profile']))
        return self.node['profile']

    @property
    def preset(self):
        if self.node['preset'] is None and self.node['instruction']['preset'] is not None:
            if self.node['instruction']['preset'] in self.configuration['study']['preset']:
                self.node['preset'] = deepcopy(self.configuration['study']['preset'][self.node['instruction']['preset']])
                if 'unique' in self.node['instruction']:
                    self.node['preset']['unique'] = self.node['instruction']['unique']
            else:
                raise ValueError('preset {} is undefined'.format(self.node['instruction']['preset']))
        return self.node['preset']

    @property
    def repertoire(self):
        if self.node['repertoire'] is None and self.preset is not None:
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
                    'head.organism name': self.preset['organism name'],
                    'head.strain': self.preset['strain'],
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
    def mongodb_query(self):
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
        return assemble(deepcopy(self.query))

    def cursor(self, collection):
        cursor = None
        if collection:
            cursor = self.resolver.database[collection].find(self.mongodb_query)

            if self.limit is not None:
                cursor.limit(self.limit)

            if self.skip is not None:
                cursor.skip(self.skip)

        return cursor

    def count(self, collection):
        return self.resolver.database[collection].count(self.mongodb_query)


class Pipeline(object):
    def __init__(self, configuration):
        self.log = logging.getLogger('Pipeline')
        self.configuration = configuration
        self.resolver = Resolver(self)
        
        for k,v in self.configuration['rss'].items():
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
        if request.action in self.configuration['interface']['implementation']:
            action = getattr(self, self.configuration['interface']['implementation'][request.action], None)
            if action is not None:
                action(request)
            else:
                raise ValueError('action {} is not implemented'.format(request.action))
        else:
            raise ValueError('action {} is not implemented'.format(request.action))

    # sample
    def sample_populate(self, request):
        count = 0
        if not request.instruction['library']:
            raise ValueError('must specify a library to populate')
        block = Block(self)
        collection = self.resolver.database['sample']
        if request.instruction['drop']:
            self.resolver.sample_drop(request.instruction['library'])
            self.resolver.analyzed_sample_drop(request.instruction['library'])

        while block.fill(request.instruction['library'], request.instruction['strain']):
            bulk = pymongo.bulk.BulkOperationBuilder(collection)
            for sample in block.buffer:
                document = sample.document
                if '_id' in document:
                    bulk.find({'_id': sample.document['_id']}).upsert().replace_one(sample.document)
                else:
                    bulk.insert(document)
            try:
                bulk.execute()
            except BulkWriteError as e:
                self.log.critical(e.details)
                raise SystemExit()
            count += block.size
            self.log.info('%s so far', count)

    def sample_drop(self, request):
        self.resolver.sample_drop(request.instruction['library'])

    def sample_view(self, request):
        request.kind = 'sample'
        cursor = request.cursor('analyzed_sample')
        for node in cursor:
            sample = Sample(self, node)
            sample.view(request.instruction['profile'])
        cursor.close()

    def sample_analyze(self, request):
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

        request.kind = 'sample'
        cursor = request.cursor('sample')
        collection = self.resolver.database['analyzed_sample']

        count = 0
        buffer = []
        for node in cursor:
            sample = Sample(self, node)
            sample.analyze()
            buffer.append(sample)
            if len(buffer) >= self.configuration['constant']['buffer size']:
                count += flush(buffer, collection)
                buffer = []
                self.log.info('%s so far', count)
        cursor.close()
        flush(buffer, collection)

    def sample_count(self, request):
        request.kind = 'sample'
        print(request.count('analyzed_sample'))

    def sample_fasta(self, request):
        request.kind = 'sample'
        cursor = request.cursor('sample')
        for node in cursor:
            sample = Sample(self, node)
            print(sample.fasta)
        cursor.close()

    def sample_fastq(self, request):
        request.kind = 'sample'
        cursor = request.cursor('sample')
        for node in cursor:
            sample = Sample(self, node)
            print(sample.fastq)
        cursor.close()

    def sample_info(self, request):
        request.kind = 'sample'
        cursor = request.cursor('analyzed_sample')
        for node in cursor:
            sample = Sample(self, node)
            sample.info(profile)
        cursor.close()

    def analyzed_sample_drop(self, request):
        self.resolver.analyzed_sample_drop(request.instruction['library'])

    # study
    def study_info(self, request):
        request.kind = 'study'
        cursor = request.cursor('study')
        for node in cursor:
            study = Study(self, node)
            study.info()
        cursor.close()

    def study_list(self, request):
        request.kind = 'study'
        cursor = request.cursor('study')
        for node in cursor:
            study = Study(self, node)
            print(str(study))
        cursor.close()

    def study_csv(self, request):
        study = self.resolver.study_resolve(request)
        if study.pivot is not None:
            study.table({ 'depth': request.instruction['depth'] })
        else:
            self.log.error('missing study cache for %s', study.name)

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
                'VH': 3,
                'DH': 2,
                'JH': 1,
            },
        }
        buffer = []

        position = 0
        for name in request.instruction['names']:
            position += 1
            order['study'][name] = position
            study = self.resolver.study_find(name)
            if study:
                study.tiled_expression(buffer)

        buffer = sorted(buffer, key=lambda x: -order['pivot'][x['pivot']])
        buffer = sorted(buffer, key=lambda x: -order['study'][x['study']])
        buffer = sorted(buffer, key=lambda x: -x['start'])

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

    # library
    def library_populate(self, request):
        for path in request.instruction['path']:
            count = 0
            with io.open(path, 'rb') as file:
                content = StringIO(file.read().decode('utf8'))
                document = json.loads(content.getvalue())
                for node in document:
                    self.resolver.library_store(node)
                    count += 1
            self.log.info('populated %d libraries', count)

    def library_info(self, request):
        request.kind = 'library'
        cursor = request.cursor('library')
        buffer = []
        for node in cursor:
            buffer.append(node)
        cursor.close()
        buffer = sorted(buffer, key=lambda x: '' if 'technical' not in x else x['technical'])
        buffer = sorted(buffer, key=lambda x: '' if 'biological' not in x else x['biological'])
        buffer = sorted(buffer, key=lambda x: '' if 'tissue' not in x else x['head']['tissue'])
        buffer = sorted(buffer, key=lambda x: '' if 'exposure' not in x else x['head']['exposure'])
        print(to_json(buffer))

    def library_list(self, request):
        title = ' | '.join (
            ['{:<9}', '{:<3}', '{:<2}', '{:<2}', '{:<30}', '{:12}', '{:8}', '{:8}', '{:<25}', '{}']
        ).format (
            'Strain', 'Exp', 'B', 'T', 'Tissue', 'Count', 'Valid', 'Prod', 'ID', 'Path'
        )
        print('| ' + title)

        separator = '-|-'.join (
            ['{:-<9}', '{:-<3}', '{:-<2}', '{:-<2}', '{:-<30}', '{:-<12}', '{:-<8}', '{:-<8}', '{:-<25}', '{:-<4}']
        ).format (
            '', '', '', '', '', '', '', '', '', ''
        )
        print('| ' + separator)

        buffer = {}
        order = []
        request.kind = 'library'
        cursor = request.cursor('library')
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
                    request = Request(self, {'library': library.id, 'profile': 'all', 'kind': 'sample'})
                    library.count = request.count('sample')

                    request = Request(self, {'library': library.id, 'profile': 'p', 'kind': 'sample'})
                    library.productive_count = request.count('analyzed_sample')

                    request = Request(self, {'library': library.id, 'profile': 'valid', 'kind': 'sample'})
                    library.valid_count = request.count('analyzed_sample')

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

        for library in order:
            print(str(library))

    # gene
    def gene_populate(self, request):
        for path in request.instruction['path']:
            count = 0
            with io.open(path, 'rb') as file:
                document = json.loads(file.read().decode('utf8'))
                for node in document:
                    self.resolver.gene_store(node)
                    count += 1
            self.log.info('populated %d genes', count)

    def gene_info(self, request):
        request.kind = 'gene'
        cursor = request.cursor('gene')
        buffer = []
        for node in cursor:
            document = node['body'].copy()
            document.update(node['head'])
            buffer.append(document)
        cursor.close()
        buffer = sorted(buffer, key=lambda x: '' if 'allele' not in x else x['allele'])
        buffer = sorted(buffer, key=lambda x: '' if 'gene' not in x else x['gene'])
        buffer = sorted(buffer, key=lambda x: '' if 'family' not in x else x['family'])
        buffer = sorted(buffer, key=lambda x: '' if 'strain' not in x else x['strain'])
        print(to_json(buffer))

    def gene_list(self, request):
        request.kind = 'gene'
        cursor = request.cursor('gene')
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

    def gene_html(self, request):
        request.kind = 'gene'
        cursor = request.cursor('gene')
        buffer = []
        for node in cursor:
            buffer.append(Gene(self, node))
        cursor.close()
        buffer = sorted(buffer, key=lambda x: '' if 'allele' not in x.body else x.body['allele'])
        buffer = sorted(buffer, key=lambda x: '' if 'gene' not in x.body else x.body['gene'])
        buffer = sorted(buffer, key=lambda x: '' if 'family' not in x.body else x.body['family'])
        buffer = sorted(buffer, key=lambda x: '' if 'strain' not in x.body else x.body['strain'])
        buffer = sorted(buffer, key=lambda x: '' if 'reference start' not in x.body else x.body['reference start'])

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

    def gene_rss(self, request):
        request.kind = 'gene'
        cursor = request.cursor('gene')
        for node in cursor:
            gene = Gene(self, node)
            gene.check_rss(request.instruction['flanking'], request.instruction['distance'])
            self.resolver.gene_save(gene)
        cursor.close()

    def gene_count(self, request):
        request.kind = 'gene'
        print(request.count('gene'))

    def gene_align(self, request):
        instruction = {
            'record': {}, 
            'total': 0, 
            'flank': request.instruction['flanking'],
            'target path': self.configuration['reference']['mouse chromosome 12']
        }
        
        request.kind = 'gene'
        cursor = request.cursor('gene')
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
                    self.log.info('gene %s aligned to multiple locations', gene.id)
            else:
                self.log.error('no satisfactory alignment found for gene %s',  gene.id)

    def gene_fasta(self, request):
        request.kind = 'gene'
        cursor = request.cursor('gene')
        buffer = []
        for node in cursor:
            buffer.append(node)
        cursor.close()
        buffer = sorted(buffer, key=lambda x: '' if 'allele' not in x else x['body']['allele'])
        buffer = sorted(buffer, key=lambda x: '' if 'gene' not in x else x['body']['gene'])
        buffer = sorted(buffer, key=lambda x: '' if 'family' not in x else x['body']['family'])
        buffer = sorted(buffer, key=lambda x: '' if 'strain' not in x else x['body']['strain'])
        for node in buffer:
            gene = Gene(self, node)
            print(gene.to_fasta(request.instruction['flanking'], self.configuration['constant']['fasta line length']))
        cursor.close()

    def gene_igblast_auxiliary(self, request):
        request.kind = 'gene'
        cursor = request.cursor('gene')
        buffer = StringIO()
        for node in cursor:
            gene = Gene(self, node)
            if gene.framed:
                buffer.write(gene.id)
                buffer.write('\t')
                buffer.write(str(gene.sequence.read_frame))
                buffer.write('\t')
                buffer.write(gene.region)
                buffer.write('\n')
        buffer.seek(0)
        print(buffer.read())

    # accession
    def accession_fasta(self, request):
        request.kind = 'accession'
        cursor = request.cursor('accession')
        for node in cursor:
            accession = Accession(self, node)
            print(accession.fasta)
        cursor.close()

    def rebuild(self, request):
        self.resolver.rebuild(request.instruction['table'])


def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    
    configuration = load_configuration()    
    command = CommandLineParser(configuration['interface'])
    if command.sectioned and command.action is None:
        command.help()

    else:
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
