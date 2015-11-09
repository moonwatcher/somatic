#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import io
import logging
import sys
import math 

iupac_nucleic_acid_notation = {
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
}

reverse_complement = {
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

class Block(object):
    def __init__(self, one_path, two_path, merged_path):
        self.log = logging.getLogger('Block')
        self.one_path = one_path
        self.two_path = two_path
        self.merged_path = merged_path
        self.fasqt_one = None
        self.fasqt_two = None
        self.fasqt_merged = None
        self.buffer = []
        self.ambiguity = self.make_ambiguity_map()
        self.ambiguity_count = {
            'nucleotide': {},
            'position': {},
            'quality': {}
        }
        self.normal = {
            'score': {},
            'phred': {}
        }
        self.score_threshold = 60

    def report(self):
        def sort_distribution(m):
            return ['{:2} {:5}'.format(p[0], p[1]) for p in sorted(m.items(), key=lambda x: x[1], reverse=True)]

        def print_normal(name):
            m = self.normal[name]
            count = 0
            mean = 0
            sd = 0

            for k,v in m.items():
                count += v
                mean += (k * v)
            if count > 0:
                mean = float(mean) / float(count)

                for k,v in m.items():
                    sd += (v * math.pow(k - mean, 2))

                sd = float(sd) / float(count)
                sd = math.sqrt(sd)
            self.log.info('\n{} distribution: count {}, mean {}, sd {}'.format(name, count, mean, sd))

        accumulative = {}
        for k,v in self.ambiguity_count['nucleotide'].items():
             for n in iupac_nucleic_acid_notation[k]['option']:
                if n not in accumulative:
                    accumulative[n] = 0
                accumulative[n] += v

        self.log.info('\nambiguity distribution\n{}'.format('\n'.join(sort_distribution(self.ambiguity_count['nucleotide']))))
        self.log.info('\nambiguity accumulative distribution\n{}'.format('\n'.join(sort_distribution(accumulative))))
        self.log.info('\nambiguity position distribution\n{}'.format('\n'.join(sort_distribution(self.ambiguity_count['position']))))
        self.log.info('\nambiguity quality distribution\n{}'.format('\n'.join(sort_distribution(self.ambiguity_count['quality']))))
        print_normal('score')
        print_normal('phred')

    def count_ambiguity(self, nucleotide, position, quality):
        if nucleotide not in self.ambiguity_count['nucleotide']:
            self.ambiguity_count['nucleotide'][nucleotide] = 0
        self.ambiguity_count['nucleotide'][nucleotide] += 1

        if position not in self.ambiguity_count['position']:
            self.ambiguity_count['position'][position] = 0
        self.ambiguity_count['position'][position] += 1

        if quality not in self.ambiguity_count['quality']:
            self.ambiguity_count['quality'][quality] = 0
        self.ambiguity_count['quality'][quality] += 1

    def count(self, read):
        if 'score' in read['picked']:
            score = read['picked']['score']
            if score not in self.normal['score']:
                self.normal['score'][score] = 0
            self.normal['score'][score] += 1

        phred = round(read['picked']['phred'])
        if phred not in self.normal['phred']:
            self.normal['phred'][phred] = 0
        self.normal['phred'][phred] += 1

    def make_ambiguity_map(self):
        solution = {}
        for i in iupac_nucleic_acid_notation.keys():
            for j in iupac_nucleic_acid_notation.keys():
                if i not in solution:
                    solution[i] = {}
                solution[i][j] = set(iupac_nucleic_acid_notation[i]['option']).union(set(iupac_nucleic_acid_notation[j]['option']))
        for i in iupac_nucleic_acid_notation.keys():
            for j in iupac_nucleic_acid_notation.keys():
                c = solution[i][j]
                for k,v in iupac_nucleic_acid_notation.items():
                    o = set(v['option'])
                    if o == c:
                        solution[i][j] = k
        return solution

    def disambiguate(self, first, second):
        return self.ambiguity[first][second]

    def open(self):
        self.fasqt_one = io.open(self.one_path, 'r')
        self.fasqt_two = io.open(self.two_path, 'r')
        self.fasqt_merged = io.open(self.merged_path, 'w')

    def close(self):
        self.fasqt_one.close()
        self.fasqt_two.close()
        self.fasqt_merged.close()

    @property
    def size(self):
        return len(self.buffer)

    @property
    def empty(self):
        return self.size == 0

    @property
    def full(self):
        return not self.size < 128

    def pick(self, read):
        if read['merged']['score'] > self.score_threshold:
            read['picked'] = read['merged']
        else:
            self.log.info(
                'not merging  {:50} S:{} O:{} Q1:{:.3} Q2:{:.3}'.format(
                    read['merged']['id'],
                    read['merged']['score'],
                    read['merged']['offset'],
                    read['one']['phred'],
                    read['two']['phred']
                )
            )
            if read['one']['phred'] > read['two']['phred']:
                read['picked'] = read['one']
            else:
                read['picked'] = self.complement(read['two'])
                for k,v in read['two'].items():
                    if k not in ['nucleotide', 'quality']:
                        read['picked'][k] = v

    def fill(self):
        self.buffer = []
        while not self.full:
            read = self.read()
            if read is not None:
                self.buffer.append(read)
                self.pick(read)
                self.count(read)
            else:
                break
        return not self.empty

    def write(self):
        for sample in self.buffer:
            self.fasqt_merged.write(sample['picked']['id'])
            self.fasqt_merged.write('\n')
            self.fasqt_merged.write(sample['picked']['nucleotide'])
            self.fasqt_merged.write('\n+\n')
            self.fasqt_merged.write(sample['picked']['quality'])
            self.fasqt_merged.write('\n')

    def read(self):
        def parse(file):
            sample = {
                'id': None,
                'nucleotide': None,
                'quality': None
            }
            state = 0
            for line in file:
                if line:
                    line = line.strip()
                    if line:
                        if state == 0 and line[0] == '@':
                            sample['id'] = line
                            state = 1
                            
                        elif state == 1:
                            sample['nucleotide'] = line
                            state = 2
                            
                        elif state == 2:
                            state = 3
                            
                        elif state == 3:
                            sample['quality'] = line
                            break
                else:
                    break
            if (sample['id'] is not None and \
                sample['nucleotide'] is not None and \
                sample['quality'] is not None and \
                len(sample['quality']) == len(sample['nucleotide'])):
                sample['length'] = len(sample['nucleotide'])
                self.average_phred(sample)
            else:
                sample = None
            return sample

        answer = { 'one': parse(self.fasqt_one), 'two': parse(self.fasqt_two) }
        if answer['one'] is not None and answer['two'] is not None:
            answer['merged'] = self.merge(answer['one'], answer['two'])
        else:
            answer = None
        return answer

    def quality(self, nucleotide):
        return ord(nucleotide) - 33

    def average_phred(self, read):
        if read and read['quality']:
            q = [ self.quality(n) for n in read['quality'] ]
            read['phred'] = float(sum(q)) / float(len(q))

    def complement(self, read):
        return {
            'nucleotide': ''.join([ reverse_complement[n] for n in read['nucleotide'] ][::-1]),
            'quality': read['quality'][::-1]
        }
    
    def align(self, one, two):
        best = None
        answer = None
        for offset in range(len(one['nucleotide'])):
            score = 0
            length = min(len(one['nucleotide']) - offset, len(two['nucleotide']))
            if length > 0:
                for i in range(length):
                    if one['nucleotide'][i + offset] == two['nucleotide'][i]:
                        score += 1
                    else:
                        score -= 1

                if best is None or score > best:
                    best = score
                    answer = offset
        return best, answer

    def merge(self, one, two):
        two = self.complement(two)
        answer = {
            'id': one['id'],
            'nucleotide': [],
            'quality': []
        }

        score_left, offset_left = self.align(one, two)
        score_right, offset_right = self.align(two, one)
        if score_left > score_right:
            offset = offset_left
            score = score_left
            start = 0
            end = max(len(one['nucleotide']), len(two['nucleotide']) + offset )
        else:
            offset = -offset_right
            score = score_right
            start = offset
            end = max(len(one['nucleotide']), len(two['nucleotide']) + offset)

        answer['offset'] = offset
        answer['score'] = score

        for position in range(start, end):
            o = None
            t = None
            if position >= 0 and position < len(one['nucleotide']):
                o = position

            if position - offset >= 0 and position - offset < len(two['nucleotide']):
                t = position - offset

            if o is None and t is not None:
                answer['nucleotide'].append(two['nucleotide'][t])
                answer['quality'].append(two['quality'][t])

            elif o is not None and t is None:
                answer['nucleotide'].append(one['nucleotide'][o])
                answer['quality'].append(one['quality'][o])

            elif o is not None and t is not None:
                if one['quality'][o] > two['quality'][t]:
                    answer['nucleotide'].append(one['nucleotide'][o])
                    answer['quality'].append(one['quality'][o])

                elif one['quality'][o] < two['quality'][t]:
                    answer['nucleotide'].append(two['nucleotide'][t])
                    answer['quality'].append(two['quality'][t])

                else:
                    nucleotide = self.disambiguate(one['nucleotide'][o], two['nucleotide'][t])
                    answer['nucleotide'].append(nucleotide)
                    answer['quality'].append(one['quality'][o])
                    if one['nucleotide'][o] != two['nucleotide'][t]:
                        self.count_ambiguity(nucleotide, position, self.quality(one['quality'][o]))
                        self.log.info(
                            'disambiguate {:58} {:3} {:2}:{}:{} --> {}'.format(
                                answer['id'],
                                position,
                                self.quality(one['quality'][o]),
                                one['nucleotide'][o],
                                two['nucleotide'][t],
                                nucleotide
                            )
                        )

        answer['nucleotide'] = ''.join(answer['nucleotide'])
        answer['quality'] = ''.join(answer['quality'])
        answer['length'] = len(answer['nucleotide'])
        self.average_phred(answer)
        return answer

    def process(self):
        self.open()
        while self.fill():
            self.write()
        self.close()

def to_json(node):
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4)

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)

    if len(sys.argv) > 3:
        one_path = sys.argv[1]
        two_path = sys.argv[2]
        merged_path = sys.argv[3]

        block = Block(one_path, two_path, merged_path)
        block.process()
        block.report()
        sys.exit(0)
    else:
        print('usage: mergeseq <first fastq> <second fastq> <merged fastq>')
        sys.exit(1)

if __name__ == '__main__':
    main()
