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

def hammig(one, two):
    distance = 0
    for i in range(len(one)):
        if one[i] != two[i]:
            distance += 1
    return distance

class SequenceDenoiser(object):
    def __init__(self, configuration=None):
        self.log = logging.getLogger('Noise')
        self.configuration = {
            'output path': None,
            'iteration': 5,
            'abundance ratio': 10,
            'read length': 145,
            'similarity threshold': 2,
            'ambiguous code': 'BDHKMNRSVWY',
        }
        self.output = None
        self.statitic = {
            'count': 0
        }
        self.unique = {}
        self.putative = []
        self.singleton = []
        self.uncertain = []
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
        self.queue.sort(key=lambda c: c['error'], reverse=False)
        self.queue.sort(key=lambda c: c['abundance'], reverse=True)

    def add(self, segment):
        if len(segment['nucleotide']) >= self.length:
            self.statitic['count'] += 1
            key = segment['nucleotide'][:self.length]
            if key in self.unique:
                cluster = self.unique[key]
                cluster['abundance'] += 1
            else:
                cluster = {
                    'key': key,
                    'id': segment['id'],
                    'abundance': 1,
                    'ambiguous': False,
                    'phred': [ ],
                    'nucleotide': set(),
                }
                self.unique[key] = cluster
                for c in self.configuration['ambiguous code']:
                    if c in key:
                        cluster['ambiguous'] = True
                        break

            cluster['phred'].append(segment['quality'])
            cluster['nucleotide'].add(segment['nucleotide'])

    def write(self):
        for key, cluster in self.unique.items():
            id = '{} {}:{}'.format(cluster['id'], cluster['abundance'], 'Y' if cluster['ambiguous'] else 'N')
            self.output.write(to_fastq(id, cluster['key'],cluster['quality']))

    def denoise(self):
        self.consolidate_quality()
        self.prepare()
        iteration = 1
        before = len(self.queue) + 1
        while iteration <= self.configuration['iteration'] and before > len(self.queue):
            self.log.info('absorption iteration %d started with %s sequences in queue', iteration, len(self.queue))
            self.bubble()
            before = len(self.queue)
            iteration += 1
        self.log.info(to_json(self.queue))

    def consolidate_quality(self):
        for cluster in self.unique.values():
            if cluster['abundance'] > 1:
                cluster['quality'] = []
                for i in range(self.length):
                    q = 0
                    for j in cluster['phred']:
                        q += phred_code_to_value(j[i])

                    cluster['quality'].append(phred_value_to_code(q / len(cluster['phred'])))
                cluster['quality'] = ''.join(cluster['quality'])
            else:
                cluster['quality'] = cluster['phred'][0][:self.length]

            cluster['error'] = expected_error(cluster['quality'])

    def prepare(self):
        self.queue = []
        self.pending = []
        for key,cluster in self.unique.items():
            if cluster['abundance'] == 1:
                cluster['state'] = 'singleton'
                self.queue.append(cluster)
            elif cluster['abundance'] >= self.configuration['abundance ratio']:
                cluster['state'] = 'putative'
                self.queue.insert(0, cluster)
            else: 
                cluster['state'] = 'pending'
                self.queue.append(cluster)
        self.sort()

    def bubble(self):
        while self.queue and self.queue[-1]['state'] != 'putative':
            cluster = self.queue.pop()
            if not self.classify(cluster):
                self.pending.append(cluster)
            else:
                self.queue.sort(key=lambda c: c['abundance'], reverse=True)

        self.queue.extend(self.pending)
        self.pending = []
        self.sort()

    def hammig(self, left, right):
        if left > right:
            if left not in self.cache:
                self.cache[left] = {}

            if right not in self.cache[left]:
                self.cache[left][right] = hammig(left, right)

            result = self.cache[left][right]
        else:
            if right not in self.cache:
                self.cache[right] = {}
            if left not in self.cache[right]:
                self.cache[right][left] = hammig(right, left)

            result = self.cache[right][left]

        return result

    def report_absorb(self, cluster, putative):
        report = []
        report.append('absorbing {} into {}'.format(cluster['id'], putative['id']))
        report.append('{:<4} '.format(putative['abundance']) + putative['key'])
        report.append('{:<4} '.format('') + putative['quality'])
        report.append('{:<4} '.format(cluster['abundance']) + cluster['key'])
        report.append('{:<4} '.format('') + cluster['quality'])
        check = []
        for i in range(len(putative['key'])):
            if putative['key'][i] == cluster['key'][i]:
                check.append('-')
            else:
                check.append('*')
        report.append('{:<4} '.format('') + ''.join(check))
        report.append('')
        sys.stderr.write('\n'.join(report))

    def absorb(self, cluster, putative):
        self.report_absorb(cluster, putative)

        # correct quality
        quality = []
        for i in range(self.length):
            q = (phred_code_to_value(putative['quality'][i]) * putative['abundance'] + phred_code_to_value(cluster['quality'][i]) * cluster['abundance']) / (putative['abundance'] + cluster['abundance'])
            quality.append(phred_value_to_code(q))
        putative['quality'] = ''.join(quality)

        putative['abundance'] += 1
        putative['phred'].append(cluster['quality'])
        putative['nucleotide'] |= cluster['nucleotide']
        if putative['abundance'] >= self.configuration['abundance ratio']:
            putative['state'] = 'putative'

    def classify(self, cluster):
        result = False
        for threshold in range(self.configuration['similarity threshold']):
            for putative in self.queue:
                distance = self.hammig(cluster['key'], putative['key'])
                ratio = putative['abundance'] / cluster['abundance']
                if ratio >= self.configuration['abundance ratio'] and distance <= threshold + 1:
                    self.absorb(cluster, putative)
                    del self.unique[cluster['key']]
                    result = True
                    break
            if result:
                break
        return result

class FastqFilter(object):
    def __init__(self, configuration=None):
        self.log = logging.getLogger('Filter')
        self.configuration = {
            'forward path': None,
            'reverse path': None,
            'merged path': None,
            'quota': 128,
            'limit': None,
            'minimum read length': 145,
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
            'retain unmerged': False,
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
        self.buffer = []
        self.output = []
        self.eof = False
        self.statitic = {
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
    def count(self):
        return self.statitic['total']

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
        self.forward = io.open(self.configuration['forward path'], 'r')
        self.reverse = io.open(self.configuration['reverse path'], 'r')
        if self.configuration['merged path'] is not None:
            self.merged = io.open(self.configuration['merged path'], 'w')
        else:
            self.merged = sys.stdout

    def close(self):
        if self.forward is not None:
            self.forward.close()

        if self.reverse is not None:
            self.reverse.close()

        if self.configuration['merged path'] is not None:
            self.merged.close()

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
            sample =  {
                'forward': self.read(self.forward),
                'reverse': self.read(self.reverse)
            }
            if sample['forward'] and sample['reverse']:
                fid = sample['forward']['id'].split()[0]
                rid = sample['reverse']['id'].split()[0]
                if fid == rid:
                    sample['id'] = sample['forward']['id'].split()[0]
                    sample['forward']['origin'] = ['forward']
                    sample['reverse']['origin'] = ['reverse']
                    self.buffer.append(sample)
                else:
                    raise ValueError('forward and reverse read don\'t match %s %s', fid, rid)

        self.process_buffer()
        self.statitic['total'] += len(self.buffer)
        return self.size > 0

    def process_buffer(self):
        for sample in self.buffer:
            self.locate_adapter(sample['forward'])
            self.locate_adapter(sample['reverse'])
            self.merge_sample(sample)
            self.pick_sample(sample)

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
            right = len(segment['nucleotide'])
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
        self.statitic['matched adapter'] += 1

        if adapter['score'] not in self.statitic['adapter score']:
            self.statitic['adapter score'][adapter['score']] = 0
        self.statitic['adapter score'][adapter['score']] += 1

        if adapter['distance'] not in self.statitic['adapter mismatch']:
            self.statitic['adapter mismatch'][adapter['distance']] = 0
        self.statitic['adapter mismatch'][adapter['distance']] += 1

        if adapter['length'] not in self.statitic['adapter overlap']:
            self.statitic['adapter overlap'][adapter['length']] = 0
        self.statitic['adapter overlap'][adapter['length']] += 1

        self.log.debug('%s found adapter %s at position %s with %s mismatch', segment['id'], adapter['nucleotide'], adapter['start'], adapter['distance'])

    def report_merge(self, sample):
        self.statitic['merged'] += 1

        score = round(sample['merged']['score'] / 10) * 10
        if score not in self.statitic['merged score']:
            self.statitic['merged score'][score] = 0
        self.statitic['merged score'][score] += 1

        mismatch = round(sample['merged']['mismatch'] / 10) * 10
        if mismatch not in self.statitic['merged mismatch']:
            self.statitic['merged mismatch'][mismatch] = 0
        self.statitic['merged mismatch'][mismatch] += 1

        overlap = round(sample['merged']['overlap'] / 10) * 10
        if overlap not in self.statitic['merged overlap']:
            self.statitic['merged overlap'][overlap] = 0
        self.statitic['merged overlap'][overlap] += 1

        offset = round(sample['merged']['offset'] / 10) * 10
        if offset not in self.statitic['merged offset']:
            self.statitic['merged offset'][offset] = 0
        self.statitic['merged offset'][offset] += 1

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
        if 'merged' in sample:
            self.trim_segment(sample['merged'])
            picked = sample['merged']

        elif self.configuration['retain unmerged']:
            self.trim_segment(sample['forward'])
            self.trim_segment(sample['reverse'])
            if sample['forward']['expected error'] < sample['forward']['expected error']:
                picked = copy.deepcopy(sample['forward'])
                self.log.info('retaining forward read %s', sample['id'])
            else:
                picked = self.complement_segment(sample['reverse'])
                self.log.info('retaining reversed read %s', sample['id'])
        else:
            self.log.info('dropping unmerged read pair %s', sample['id'])

        # drop reads that are shorter than minimum
        if picked is not None and picked['length'] < self.configuration['minimum read length']:
            self.log.info('dropping short read %s of length %d', picked['id'], picked['length'])
            picked = None

        if picked is not None:
            picked['id'] = sample['id']
            self.output.append(picked)

    def trim_segment(self, segment):
        nucleotide = segment['nucleotide']
        quality = segment['quality']
        length = len(nucleotide)
        window = self.configuration['trimming window']
        threshold = self.configuration['trimming error threshold']

        if self.configuration['trim end']:
            while expected_error(quality[-window:]) > threshold and length >= self.configuration['minimum read length']:
                nucleotide = nucleotide[:-1]
                quality = quality[:-1]
                length -= 1

        if self.configuration['trim start']:
            while expected_error(quality[:window]) > threshold and length >= self.configuration['minimum read length']:
                nucleotide = nucleotide[1:]
                quality = quality[1:]
                length -= 1

        segment['nucleotide'] = nucleotide
        segment['quality'] = quality
        segment['length'] = length
        segment['expected error'] = expected_error(segment['quality'])
        for term in [
            'effective length',
            'effective start',
            'effective end'
        ]:
            if term in segment:
                del segment[term]

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
    fastqFilter = FastqFilter(
        {
            'forward path': configuration['forward path'],
            'reverse path': configuration['reverse path'],
            'limit': configuration['limit'],
        }
    )
    sequenceDenoiser = SequenceDenoiser(
        {
            'output path': configuration['output path'],
        }
    )

    sequenceDenoiser.open()

    fastqFilter.open()
    while fastqFilter.fill() and not fastqFilter.enough:
        # fastqFilter.write()
        for sample in fastqFilter.output:
            sequenceDenoiser.add(sample)
    fastqFilter.close()

    sequenceDenoiser.denoise()
    sequenceDenoiser.write()
    sequenceDenoiser.close()

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    if len(sys.argv) > 2:
        configuration = {
            'forward path': sys.argv[1],
            'reverse path': sys.argv[2],
            'output path': None,
            'limit': None,
        }
        if len(sys.argv) > 3:
            configuration['output path'] = sys.argv[3]

        process(configuration)
        sys.exit(0)
    else:
        print('usage: mergeseq <first fastq> <second fastq> <merged fastq>')
        sys.exit(1)

if __name__ == '__main__':
    main()
