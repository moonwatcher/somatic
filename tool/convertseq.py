#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Extract a sequence by coordinates from a FASTA file
# Author: Lior Galanti < lior.galanti@nyu.edu >
# NYU Center for Genetics and System Biology 2015
#
# crop is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation;
# either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.

import sys
import json
import logging
from argparse import ArgumentParser

configuration = {
    'interface': {
        'global': {
            'argument': [
                'from',
                'to',
                'version', 
                'verbosity'
            ]
        }, 
        'instruction': {
            'description': 'Lior Galanti lior.galanti@nyu.edu NYU Center for Genomics & Systems Biology'
        }, 
        'prototype': {
            'from': {
                'flag': [
                    '-f', 
                    '--from'
                ], 
                'parameter': {
                    'choices': [
                        'fasta', 
                        'fastq'
                    ], 
                    'default': 'fastq', 
                    'dest': 'from', 
                    'help': 'input format'
                }
            }, 
            'to': {
                'flag': [
                    '-t', 
                    '--to'
                ], 
                'parameter': {
                    'choices': [
                        'fasta', 
                        'fastq'
                    ], 
                    'default': 'fasta', 
                    'dest': 'to', 
                    'help': 'output format'
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
                    'version': '%[prog]s 1.0'
                }
            }
        }
    }, 
    'iupac amino acid notation': {
        '*': {
            'code': 'stop', 
            'name': 'Stop codon'
        }, 
        'A': {
            'charge': 1.8, 
            'class code': 'N', 
            'class name': 'nonpolar', 
            'code': 'ala', 
            'name': 'Alanine', 
            'weight': 89
        }, 
        'B': {
            'code': 'asx', 
            'name': 'Aspartic acid or Asparagine'
        }, 
        'C': {
            'charge': 2.5, 
            'class code': 'P', 
            'class name': 'polar', 
            'code': 'cys', 
            'name': 'Cysteine', 
            'weight': 121
        }, 
        'D': {
            'charge': -3.5, 
            'class code': 'A', 
            'class name': 'acidic', 
            'code': 'asp', 
            'name': 'Aspartic acid', 
            'weight': 133
        }, 
        'E': {
            'charge': -3.5, 
            'class code': 'A', 
            'class name': 'acidic', 
            'code': 'glu', 
            'name': 'Glutamic acid', 
            'weight': 147
        }, 
        'F': {
            'charge': 2.8, 
            'class code': 'N', 
            'class name': 'nonpolar', 
            'code': 'phe', 
            'name': 'Phenylalanine', 
            'weight': 165
        }, 
        'G': {
            'charge': -0.4, 
            'class code': 'N', 
            'class name': 'nonpolar', 
            'code': 'gly', 
            'name': 'Glycine', 
            'weight': 75
        }, 
        'H': {
            'charge': -3.2, 
            'class code': 'B', 
            'class name': 'basic', 
            'code': 'his', 
            'name': 'Histidine', 
            'weight': 155
        }, 
        'I': {
            'charge': 4.5, 
            'class code': 'N', 
            'class name': 'nonpolar', 
            'code': 'ile', 
            'name': 'Isoleucine', 
            'weight': 131
        }, 
        'J': {
            'code': 'xle', 
            'name': 'Leucine or Isoleucine'
        }, 
        'K': {
            'charge': -3.9, 
            'class code': 'B', 
            'class name': 'basic', 
            'code': 'lys', 
            'name': 'Lysine', 
            'weight': 146
        }, 
        'L': {
            'charge': 3.8, 
            'class code': 'N', 
            'class name': 'nonpolar', 
            'code': 'leu', 
            'name': 'Leucine', 
            'weight': 131
        }, 
        'M': {
            'charge': 1.9, 
            'class code': 'N', 
            'class name': 'nonpolar', 
            'code': 'met', 
            'name': 'Methionine', 
            'weight': 149
        }, 
        'N': {
            'charge': -3.5, 
            'class code': 'P', 
            'class name': 'polar', 
            'code': 'asn', 
            'name': 'Asparagine', 
            'weight': 132
        }, 
        'P': {
            'charge': -1.6, 
            'class code': 'N', 
            'class name': 'nonpolar', 
            'code': 'pro', 
            'name': 'Proline', 
            'weight': 115
        }, 
        'Q': {
            'charge': -3.5, 
            'class code': 'P', 
            'class name': 'polar', 
            'code': 'gln', 
            'name': 'Glutamine', 
            'weight': 146
        }, 
        'R': {
            'charge': -4.5, 
            'class code': 'B', 
            'class name': 'basic', 
            'code': 'arg', 
            'name': 'Arginine', 
            'weight': 174
        }, 
        'S': {
            'charge': -0.8, 
            'class code': 'P', 
            'class name': 'polar', 
            'code': 'ser', 
            'name': 'Serine', 
            'weight': 105
        }, 
        'T': {
            'charge': -0.7, 
            'class code': 'P', 
            'class name': 'polar', 
            'code': 'thr', 
            'name': 'Threonine', 
            'weight': 119
        }, 
        'V': {
            'charge': 4.2, 
            'class code': 'N', 
            'class name': 'nonpolar', 
            'code': 'val', 
            'name': 'Valine', 
            'weight': 117
        }, 
        'W': {
            'charge': -0.9, 
            'class code': 'N', 
            'class name': 'nonpolar', 
            'code': 'trp', 
            'name': 'Tryptophan', 
            'weight': 204
        }, 
        'X': {
            'code': 'xaa', 
            'name': 'Any amino acid'
        }, 
        'Y': {
            'charge': -1.3, 
            'class code': 'P', 
            'class name': 'polar', 
            'code': 'tyr', 
            'name': 'Tyrosine', 
            'weight': 181
        }, 
        'Z': {
            'code': 'glx', 
            'name': 'Glutamine or Glutamic acid'
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
    'nucleic to amino': {
        'AAA': 'K', 
        'AAC': 'N', 
        'AAG': 'K', 
        'AAT': 'N', 
        'ACA': 'T', 
        'ACC': 'T', 
        'ACG': 'T', 
        'ACT': 'T', 
        'AGA': 'R', 
        'AGC': 'S', 
        'AGG': 'R', 
        'AGT': 'S', 
        'ATA': 'I', 
        'ATC': 'I', 
        'ATG': 'M', 
        'ATT': 'I', 
        'CAA': 'Q', 
        'CAC': 'H', 
        'CAG': 'Q', 
        'CAT': 'H', 
        'CCA': 'P', 
        'CCC': 'P', 
        'CCG': 'P', 
        'CCT': 'P', 
        'CGA': 'R', 
        'CGC': 'R', 
        'CGG': 'R', 
        'CGT': 'R', 
        'CTA': 'L', 
        'CTC': 'L', 
        'CTG': 'L', 
        'CTT': 'L', 
        'GAA': 'E', 
        'GAC': 'D', 
        'GAG': 'E', 
        'GAT': 'D', 
        'GCA': 'A', 
        'GCC': 'A', 
        'GCG': 'A', 
        'GCT': 'A', 
        'GGA': 'G', 
        'GGC': 'G', 
        'GGG': 'G', 
        'GGT': 'G', 
        'GTA': 'V', 
        'GTC': 'V', 
        'GTG': 'V', 
        'GTT': 'V', 
        'TAA': '*', 
        'TAC': 'Y', 
        'TAG': '*', 
        'TAT': 'Y', 
        'TCA': 'S', 
        'TCC': 'S', 
        'TCG': 'S', 
        'TCT': 'S', 
        'TGA': '*', 
        'TGC': 'C', 
        'TGG': 'W', 
        'TGT': 'C', 
        'TTA': 'L', 
        'TTC': 'F', 
        'TTG': 'L', 
        'TTT': 'F'
    }, 
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
}

log_levels = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

def to_json(node):
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4)

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

class Pipeline(object):
    def __init__(self, node):
        self.node = node
        self.collection = None
        self.line_length = 120
        self.capacity = 128
        self.sample = None

    @property
    def size(self):
        return len(self.collection)

    @property
    def full(self):
        return self.size >= self.capacity

    @property
    def empty(self):
        return self.size == 0

    @property
    def fastq(self):
        buffer = []
        for record in self.collection:
            buffer.append('@{}\n{}\n+\n{}'.format(record['D'], record['D'], record['Q']))
        return '\n'.join(buffer)

    @property
    def fasta(self):
        buffer = []
        for record in self.collection:
            buffer.append('>{}'.format(record['D']))
            length = len(record['S'])
            begin = 0
            end = 0
            while end < length:
                end = min(begin + self.line_length, length)
                buffer.append(record['S'][begin:end])
                begin = end
        return '\n'.join(buffer)

    def execute(self):
        if self.node['from'] == 'fastq' and self.node['to'] == 'fasta':
            self.fastq_to_fasta()
        elif self.node['from'] == 'fasta' and self.node['to'] == 'fastq':
            self.fasta_to_fastq()

    def read_fastq(self):
        self.collection = []
        state = 0
        for line in sys.stdin:
            if not line:
                break
            else:
                line = line.strip()
                if line:
                    if state == 0 and line[0] == '@':
                        self.sample = { 'D': line[1:] }
                        state = 1

                    elif state == 1:
                        if self.sample is not None:
                            self.sample['S'] = line
                        state = 2

                    elif state == 2:
                        state = 3

                    elif state == 3:
                        if self.sample is not None:
                            self.sample['Q'] = line
                            self.collection.append(self.sample)
                            self.sample = None
                            state = 0
                            if self.full:
                                break
        return not self.empty

    def read_fasta(self):
        self.collection = []
        for line in sys.stdin:
            if not line:
                break
            else:
                line = line.strip()
                if line:
                    if line[0] == '>':
                        if self.sample is not None:
                            self.sample['S'] = ''.join(self.sample['S'])
                            self.collection.append(self.sample)

                        self.sample = { 'D': line[1:], 'S': [] }
                        if self.full:
                            break

                    elif sample is not None:
                        self.sample['S'].append(line)

        return not self.empty

    def fasta_to_fastq(self):
        while self.read_fasta():
            print(self.fastq)

    def fastq_to_fasta(self):
        while self.read_fastq():
            print(self.fasta)

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    command = CommandLineParser(configuration['interface'])
    logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])
    pipeline = Pipeline(command.instruction)
    try:
        pipeline.execute()

    except ValueError as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)

    except(KeyboardInterrupt, SystemExit) as e:
        pipeline.close()
        sys.exit(1)

    sys.exit(0)

if __name__ == '__main__':
    main()