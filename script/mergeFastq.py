#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Quality filtering, merging and light dereplication for IGH sequencing
# Author: Lior Galanti < lior.galanti@nyu.edu >
# NYU Center for Genetics and System Biology 2015
#
# complement is free software; you can redistribute it and/or modify it under the terms of
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
# IUPAC substitution matrix ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4

import json
import io
import logging
import sys
import math 
import copy
import hashlib

def to_json(node):
    def handler(o):
        if isinstance(o, set):
            return list(o)
        else:
            return o
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4, default=handler)

def phred_code_to_value(code):
    return ord(code) - 33

def phred_value_to_code(value):
    return chr(round(value) + 33)

def phred_to_probability(score):
    return pow(10, -(ord(score) - 33) / 10)

def probability_to_phred(probability):
    return chr(round(-10 * math.log(probability, 10)) + 33)

def expected_error(quality):
    return sum([ phred_to_probability(q) for q in quality ])

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

def to_fastq(id, nucleotide, quality):
    return '{}\n{}\n+\n{}\n'.format(id, nucleotide, quality)

def segment_to_fastq(segment):
    return '{}\n{}\n+\n{}\n'.format(segment['id'], segment['nucleotide'], segment['quality'])

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

class SequenceDenoiser(object):
    def __init__(self, configuration=None):
        self.log = logging.getLogger('Noise')
        self.configuration = {
            'output path': None,
            'iteration': 5,
            'abundance ratio': 10,
            'read length': 133,
            'similarity threshold': 2,
            'ambiguous code': 'BDHKMNRSVWY',
        }
        self.output = None
        self.statistic = {
            'count': 0
        }
        self.cache = {}
        if configuration is not None:
            for k,v in configuration.items():
                self.configuration[k] = v

    def open(self):
        if self.configuration['output path'] is not None:
            self.output = io.open(self.configuration['output path'], 'w')
        else:
            self.output = sys.stdout

    def close(self):
        if self.configuration['output path'] is not None:
            self.output.close()

    @property
    def length(self):
        return self.configuration['read length']

    def sort(self):
        self.log.info('sorting %s reads in queue', len(self.queue))
        self.queue.sort(key=lambda c: c['error'], reverse=False)
        self.queue.sort(key=lambda c: c['abundance'], reverse=True)

    def add(self, segment):
        if len(segment['nucleotide']) >= self.length:
            self.statistic['count'] += 1
            key = segment['nucleotide'][:self.length]

            id, comment = segment['id'].split(' ')
            abundance, ambiguous = comment.split(':')
            if abundance is not None:
                try:
                    abundance = int(abundance)
                except ValueError:
                    abundance = 1
            else:
                abundance = 1

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
                for c in self.configuration['ambiguous code']:
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

    def denoise(self):
        self.prepare()
        iteration = 1
        before = len(self.queue) + 1
        while iteration <= self.configuration['iteration'] and before > len(self.queue):
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
                    q = sum([phred_code_to_value(j[i]) for j in cluster['quality']])
                    quality.append(phred_value_to_code(q / count))
                cluster['quality'] = ''.join(quality)
            else:
                cluster['quality'] = cluster['quality'][0][:self.length]
            cluster['error'] = expected_error(cluster['quality'])

            if cluster['abundance'] == 1:
                cluster['state'] = 'singleton'
            elif cluster['abundance'] >= self.configuration['abundance ratio']:
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
            quality.append(phred_value_to_code((phred_code_to_value(p) * pa + phred_code_to_value(c) * ca) / (pa + ca)))
        putative['quality'] = ''.join(quality)

        putative['abundance'] += 1
        putative['nucleotide'] |= cluster['nucleotide']
        if putative['abundance'] >= self.configuration['abundance ratio']:
            putative['state'] = 'putative'

    def classify(self, cluster):
        result = False
        factor = cluster['abundance'] * self.configuration['abundance ratio']
        for threshold in range(self.configuration['similarity threshold']):
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

class FastqFilteringReader(object):
    def __init__(self, configuration=None):
        self.log = logging.getLogger('Filter')
        self.configuration = {
            'forward path': None,
            'reverse path': None,
            'paired': False,
            'filtering': True,
            'quota': 128,
            'limit': None,
            'minimum read length': 133,
            'trim start': False,
            'trim end': True,
            'trimming window': 6,
            'trimming error threshold': 1,
            'adapter mismatch': 3,
            'adapter overlap': 6,
            'adapter score threshold': 40,
            'reverse adapter scan': True,
            'adapter': [
                'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                'TACACTCTTTCCCTACACGACGCTCTTCCGATCT',
            ],
            'retain short': True,
            'retain unmerged': True,
            'error threshold': 2,
            'merge overlap threshold': 70,
            'merge score threshold': 400,
            'substitution quality correction': 28,
            'nucleotide substitution': {
                'A': {
                    'A': 5,
                    'B': -4,
                    'C': -4,
                    'D': -1,
                    'G': -4,
                    'H': -1,
                    'K': -4,
                    'M': 1,
                    'N': -2,
                    'R': 1,
                    'S': -4,
                    'T': -4,
                    'V': -1,
                    'W': 1,
                    'Y': -4
                },
                'B': {
                    'A': -4,
                    'B': -1,
                    'C': -1,
                    'D': -2,
                    'G': -1,
                    'H': -2,
                    'K': -1,
                    'M': -3,
                    'N': -1,
                    'R': -3,
                    'S': -1,
                    'T': -1,
                    'V': -2,
                    'W': -3,
                    'Y': -1
                },
                'C': {
                    'A': -4,
                    'B': -1,
                    'C': 5,
                    'D': -4,
                    'G': -4,
                    'H': -1,
                    'K': -4,
                    'M': 1,
                    'N': -2,
                    'R': -4,
                    'S': 1,
                    'T': -4,
                    'V': -1,
                    'W': -4,
                    'Y': 1
                },
                'D': {
                    'A': -1,
                    'B': -2,
                    'C': -4,
                    'D': -1,
                    'G': -1,
                    'H': -2,
                    'K': -1,
                    'M': -3,
                    'N': -1,
                    'R': -1,
                    'S': -3,
                    'T': -1,
                    'V': -2,
                    'W': -1,
                    'Y': -3
                },
                'G': {
                    'A': -4,
                    'B': -1,
                    'C': -4,
                    'D': -1,
                    'G': 5,
                    'H': -4,
                    'K': 1,
                    'M': -4,
                    'N': -2,
                    'R': 1,
                    'S': 1,
                    'T': -4,
                    'V': -1,
                    'W': -4,
                    'Y': -4
                },
                'H': {
                    'A': -1,
                    'B': -2,
                    'C': -1,
                    'D': -2,
                    'G': -4,
                    'H': -1,
                    'K': -3,
                    'M': -1,
                    'N': -1,
                    'R': -3,
                    'S': -3,
                    'T': -1,
                    'V': -2,
                    'W': -1,
                    'Y': -1
                },
                'K': {
                    'A': -4,
                    'B': -1,
                    'C': -4,
                    'D': -1,
                    'G': 1,
                    'H': -3,
                    'K': -1,
                    'M': -4,
                    'N': -1,
                    'R': -2,
                    'S': -2,
                    'T': 1,
                    'V': -3,
                    'W': -2,
                    'Y': -2
                },
                'M': {
                    'A': 1,
                    'B': -3,
                    'C': 1,
                    'D': -3,
                    'G': -4,
                    'H': -1,
                    'K': -4,
                    'M': -1,
                    'N': -1,
                    'R': -2,
                    'S': -2,
                    'T': -4,
                    'V': -1,
                    'W': -2,
                    'Y': -2
                },
                'N': {
                    'A': -2,
                    'B': -1,
                    'C': -2,
                    'D': -1,
                    'G': -2,
                    'H': -1,
                    'K': -1,
                    'M': -1,
                    'N': -1,
                    'R': -1,
                    'S': -1,
                    'T': -2,
                    'V': -1,
                    'W': -1,
                    'Y': -1
                },
                'R': {
                    'A': 1,
                    'B': -3,
                    'C': -4,
                    'D': -1,
                    'G': 1,
                    'H': -3,
                    'K': -2,
                    'M': -2,
                    'N': -1,
                    'R': -1,
                    'S': -2,
                    'T': -4,
                    'V': -1,
                    'W': -2,
                    'Y': -4
                },
                'S': {
                    'A': -4,
                    'B': -1,
                    'C': 1,
                    'D': -3,
                    'G': 1,
                    'H': -3,
                    'K': -2,
                    'M': -2,
                    'N': -1,
                    'R': -2,
                    'S': -1,
                    'T': -4,
                    'V': -1,
                    'W': -4,
                    'Y': -2
                },
                'T': {
                    'A': -4,
                    'B': -1,
                    'C': -4,
                    'D': -1,
                    'G': -4,
                    'H': -1,
                    'K': 1,
                    'M': -4,
                    'N': -2,
                    'R': -4,
                    'S': -4,
                    'T': 5,
                    'V': -4,
                    'W': 1,
                    'Y': 1
                },
                'V': {
                    'A': -1,
                    'B': -2,
                    'C': -1,
                    'D': -2,
                    'G': -1,
                    'H': -2,
                    'K': -3,
                    'M': -1,
                    'N': -1,
                    'R': -1,
                    'S': -1,
                    'T': -4,
                    'V': -1,
                    'W': -3,
                    'Y': -3
                },
                'W': {
                    'A': 1,
                    'B': -3,
                    'C': -4,
                    'D': -1,
                    'G': -4,
                    'H': -1,
                    'K': -2,
                    'M': -2,
                    'N': -1,
                    'R': -2,
                    'S': -4,
                    'T': 1,
                    'V': -3,
                    'W': -1,
                    'Y': -2
                },
                'Y': {
                    'A': -4,
                    'B': -1,
                    'C': 1,
                    'D': -3,
                    'G': -4,
                    'H': -1,
                    'K': -2,
                    'M': -2,
                    'N': -1,
                    'R': -4,
                    'S': -2,
                    'T': 1,
                    'V': -3,
                    'W': -2,
                    'Y': -1
                }
            },
            'iupac nucleic acid notation': {
                'A': {
                    'name': 'Adenine', 
                    'option': [
                        'A'
                    ], 
                    'reverse': 'T'
                }, 
                'B': {
                    'name': 'Not A', 
                    'option': [
                        'G', 
                        'T', 
                        'C'
                    ], 
                    'reverse': 'V'
                }, 
                'C': {
                    'name': 'Cytosine', 
                    'option': [
                        'C'
                    ], 
                    'reverse': 'G'
                }, 
                'D': {
                    'name': 'Not C', 
                    'option': [
                        'G', 
                        'A', 
                        'T'
                    ], 
                    'reverse': 'H'
                }, 
                'G': {
                    'name': 'Guanine', 
                    'option': [
                        'G'
                    ], 
                    'reverse': 'C'
                }, 
                'H': {
                    'name': 'Not G', 
                    'option': [
                        'A', 
                        'C', 
                        'T'
                    ], 
                    'reverse': 'D'
                }, 
                'K': {
                    'name': 'Keto', 
                    'option': [
                        'G', 
                        'T'
                    ], 
                    'reverse': 'M'
                }, 
                'M': {
                    'name': 'Amino', 
                    'option': [
                        'A', 
                        'C'
                    ], 
                    'reverse': 'K'
                }, 
                'N': {
                    'name': 'Any', 
                    'option': [
                        'A', 
                        'G', 
                        'C', 
                        'T'
                    ], 
                    'reverse': 'N'
                }, 
                'R': {
                    'name': 'Purine', 
                    'option': [
                        'G', 
                        'A'
                    ], 
                    'reverse': 'Y'
                }, 
                'S': {
                    'name': 'Strong bonds', 
                    'option': [
                        'G', 
                        'C'
                    ], 
                    'reverse': 'S'
                }, 
                'T': {
                    'name': 'Thymine', 
                    'option': [
                        'T'
                    ], 
                    'reverse': 'A'
                }, 
                'V': {
                    'name': 'Not T', 
                    'option': [
                        'G', 
                        'C', 
                        'A'
                    ], 
                    'reverse': 'B'
                }, 
                'W': {
                    'name': 'Weak bonds', 
                    'option': [
                        'A', 
                        'T'
                    ], 
                    'reverse': 'W'
                }, 
                'Y': {
                    'name': 'Pyrimidine', 
                    'option': [
                        'T', 
                        'C'
                    ], 
                    'reverse': 'R'
                }
            },
            'reverse complement': {
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
        }
        self.configuration['substitution quality correction'] = float(self.configuration['substitution quality correction']) * 2.0
        self.forward = None
        self.reverse = None
        self.merged = None
        self.buffer = None
        self.output = None
        self.eof = False
        self.statistic = {
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
        if configuration is not None:
            for k,v in configuration.items():
                self.configuration[k] = v
        self.load()

    def load(self):
        # make sure all reverse complements of adapters are present
        adapter = set()
        for sequence in self.configuration['adapter']:
            adapter.add(sequence)
            adapter.add(self.complement_sequence(sequence))
        self.configuration['adapter'] = list(adapter)

        self.configuration['ambiguity'] = {}
        for i in self.configuration['iupac nucleic acid notation'].keys():
            for j in self.configuration['iupac nucleic acid notation'].keys():
                if i not in self.configuration['ambiguity']:
                    self.configuration['ambiguity'][i] = {}
                self.configuration['ambiguity'][i][j] = set(self.configuration['iupac nucleic acid notation'][i]['option']).union(set(self.configuration['iupac nucleic acid notation'][j]['option']))
        for i in self.configuration['iupac nucleic acid notation'].keys():
            for j in self.configuration['iupac nucleic acid notation'].keys():
                c = self.configuration['ambiguity'][i][j]
                for k,v in self.configuration['iupac nucleic acid notation'].items():
                    o = set(v['option'])
                    if o == c:
                        self.configuration['ambiguity'][i][j] = k

    @property
    def limit(self):
        return self.configuration['limit']

    @property
    def quota(self):
        return self.configuration['quota']

    @property
    def paired(self):
        return self.configuration['paired']

    @property
    def filtering(self):
        return self.configuration['filtering']

    @property
    def count(self):
        return self.statistic['total']

    @property
    def size(self):
        return len(self.buffer)

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
        return self.configuration['adapter']

    def open(self):
        self.buffer = []
        self.output = []
        if self.configuration['forward path'] is not None:
            self.forward = io.open(self.configuration['forward path'], 'r')
        else:
            self.forward = sys.stdin

        if self.configuration['reverse path'] is not None:
            self.reverse = io.open(self.configuration['reverse path'], 'r')

        if self.forward is not None and self.reverse is not None:
            self.configuration['paired'] = True

        self.merged = sys.stdout

    def close(self):
        if self.configuration['forward path'] is not None:
            self.forward.close()

        if self.configuration['reverse path'] is not None:
            self.reverse.close()

        self.buffer = None
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
                        segment = { 'id': line }
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
            self.eof = True

        return segment

    def write(self):
        for segment in self.output:
            self.merged.write(to_fastq(segment['id'], segment['nucleotide'],segment['quality']))

    def fill(self):
        self.buffer = []
        self.output = []
        while not self.full and not self.eof:
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
                self.buffer.append(sample)

        self.process_buffer()
        self.statistic['total'] += len(self.buffer)
        return self.size > 0

    def process_buffer(self):
        if self.filtering:
            for sample in self.buffer:
                self.locate_adapter(sample['forward'])
                if self.paired:
                    self.locate_adapter(sample['reverse'])
                    self.merge_sample(sample)
                self.pick_sample(sample)
        else:
            for sample in self.buffer:
                self.place_in_output(sample['forward'])

    def locate_adapter(self, segment):
        taken = False
        for index, adapter in enumerate(self.adapter):
            if self.find_adapter(segment, adapter):
                self.promote_adapter(index)
                taken = True
                break

    def find_adapter(self, segment, adapter):
        taken = False
        if self.configuration['reverse adapter scan']:
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
        while end - left >= self.configuration['adapter overlap']:
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
            if length >= self.configuration['adapter overlap']:
                for i in range(length):
                    score += self.score_base_alignment(adapter[i], sequence[i], quality[i], quality[i])
                    if adapter[i] != sequence[i]:
                        distance += 1
                        if distance > self.configuration['adapter mismatch']:
                            break

                if  distance <= self.configuration['adapter mismatch'] and \
                    score >= self.configuration['adapter score threshold']:
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

        if picked is not None and picked['score'] >= self.configuration['merge score threshold']:
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
        picked = None
        if self.paired:
            if 'merged' in sample:
                self.trim_segment(sample['merged'])
                picked = sample['merged']

            elif self.configuration['retain unmerged']:
                self.trim_segment(sample['forward'])
                self.trim_segment(sample['reverse'])
                if  sample['forward']['expected error'] < sample['reverse']['expected error'] and \
                    sample['forward']['length'] >= self.configuration['minimum read length']:
                    picked = sample['forward']

                elif sample['reverse']['length'] >= self.configuration['minimum read length']:
                    picked = self.complement_segment(sample['reverse'])
                    # if we use the reverse read we need to directionally trim it again
                    self.trim_segment(picked)
            else:
                self.log.info('dropping unmerged read pair {}\n{}\n{}'.format(sample['id'], segment_to_fastq(sample['forward']), segment_to_fastq(sample['reverse'])))
        else:
            self.trim_segment(sample['forward'])
            picked = copy.deepcopy(sample['forward'])

        # drop reads that are shorter than minimum
        if picked is not None and picked['length'] < self.configuration['minimum read length']:
            self.log.info('dropping short read {} of length {}\n{}'.format(picked['id'], picked['length'], segment_to_fastq(picked)))
            picked = None

        if picked is not None:
            effective_error = expected_error(picked['quality'][:self.configuration['minimum read length']])
            if effective_error > self.configuration['error threshold']:
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
        width = self.configuration['trimming window']
        threshold = self.configuration['trimming error threshold']

        if self.configuration['trim end']:
            window = segment['quality'][-width:]
            position = 0
            while expected_error(window) > threshold and length >= self.configuration['minimum read length']:
                length -= 1
                position -= 1
                window = segment['quality'][-width + position:position]
            segment['quality'] = segment['quality'][-position:]
            segment['nucleotide'] = segment['nucleotide'][-position:]

        if self.configuration['trim start']:
            window = segment['quality'][:width]
            position = 0
            while expected_error(window) > threshold and length >= self.configuration['minimum read length']:
                length -= 1
                position += 1
                window = segment['quality'][position:width + position]
            segment['quality'] = segment['quality'][position:]
            segment['nucleotide'] = segment['nucleotide'][position:]

        segment['length'] = len(segment['nucleotide'])
        segment['expected error'] = expected_error(segment['quality'])

    def score_base_alignment(self, ni, nj, qi, qj):
        score = self.configuration['nucleotide substitution'][ni][nj]

        # if qualities were supplied correct the score
        if qi is not None and qj is not None:
            correction =  (phred_code_to_value(qi) + phred_code_to_value(qj)) / self.configuration['substitution quality correction']
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
            if length < self.configuration['merge overlap threshold']:
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
        return self.configuration['ambiguity'][i][j]

    def complement_segment(self, segment):
        reversed = {
            'id': segment['id'],
            'nucleotide': ''.join([ self.configuration['reverse complement'][n] for n in segment['nucleotide'] ][::-1]),
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

    def complement_sequence(self, sequence):
        return ''.join([ self.configuration['reverse complement'][n] for n in sequence ][::-1])

def process(configuration):
    fastqFilter = FastqFilteringReader(configuration)
    sequenceDenoiser = SequenceDenoiser(configuration)

    sequenceDenoiser.open()
    fastqFilter.open()

    while fastqFilter.fill() and not fastqFilter.enough:
        # fastqFilter.write()
        for sample in fastqFilter.output:
            sequenceDenoiser.add(sample)

    fastqFilter.close()
    sys.stderr.write(to_json(fastqFilter.statistic))
    fastqFilter = None

    sequenceDenoiser.denoise()
    sequenceDenoiser.write()
    sequenceDenoiser.close()

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    configuration = {
        'limit': None,
        'filtering': False,
    }

    if len(sys.argv) > 1:
        configuration['forward path'] = sys.argv[1]

    if len(sys.argv) > 2:
        configuration['reverse path'] = sys.argv[2],

    process(configuration)
    sys.exit(0)

if __name__ == '__main__':
    main()
