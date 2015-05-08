#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Reverse complement an IUPAC encoded DNA string
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

import sys

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

for line in sys.stdin:
    
    try:
        print(''.join([complement[b] for b in list(line.strip())][::-1]))
        
    except KeyError as e:
        print(
            '{} is not a valid IUAPAC nucleotide notation.'.format(e) +
            'see http://en.wikipedia.org/wiki/Nucleic_acid_notation'
        )
