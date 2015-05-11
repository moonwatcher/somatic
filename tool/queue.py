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
import json

LOG = '/volume/albireo/canal/thesis/populate'

def to_json(node):
    def handler(o):
        result = None
        if isinstance(o, datetime):
            result = o.isoformat()
        if isinstance(o, ObjectId):
            result = str(o)
        if isinstance(o, set):
            result = list(o)
        return result

    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4, default=handler)


libraries = [
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_fo_1_a.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 1,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_fo_1_b.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 1,
        'technical repetition': 2,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_fo_1_c.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 1,
        'technical repetition': 3,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_fo_2_a.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 2,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_fo_2_b.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 2,
        'technical repetition': 2,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_fo_3_a.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 3,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_fo_3_b.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 3,
        'technical repetition': 2,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_fo_3_c.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 3,
        'technical repetition': 3,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_fo_4_a.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 4,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_fo_4_b.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 4,
        'technical repetition': 2,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_fo_5_a.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 5,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_fo_5_b.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 5,
        'technical repetition': 2,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_fo_5_c.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 5,
        'technical repetition': 3,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_fo_5_d.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 5,
        'technical repetition': 4,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_b1a_1.fastq.gz',
        'tissue': 'peritoneal cavity b1a',
        'biological repetition': 1,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_b1a_2.fastq.gz',
        'tissue': 'peritoneal cavity b1a',
        'biological repetition': 2,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_b1a_3.fastq.gz',
        'tissue': 'peritoneal cavity b1a',
        'biological repetition': 3,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_b1a_4.fastq.gz',
        'tissue': 'peritoneal cavity b1a',
        'biological repetition': 4,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_b1a_5.fastq.gz',
        'tissue': 'peritoneal cavity b1a',
        'biological repetition': 5,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_mz_1.fastq.gz',
        'tissue': 'marginal zone spleen',
        'biological repetition': 1,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_mz_2.fastq.gz',
        'tissue': 'marginal zone spleen',
        'biological repetition': 2,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_mz_3.fastq.gz',
        'tissue': 'marginal zone spleen',
        'biological repetition': 3,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_mz_4.fastq.gz',
        'tissue': 'marginal zone spleen',
        'biological repetition': 4,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_mz_5.fastq.gz',
        'tissue': 'marginal zone spleen',
        'biological repetition': 5,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_preb_1.fastq.gz',
        'tissue': 'bone marrow pre b',
        'biological repetition': 1,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_preb_2.fastq.gz',
        'tissue': 'bone marrow pre b',
        'biological repetition': 2,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_conv_preb_3.fastq.gz',
        'tissue': 'bone marrow pre b',
        'biological repetition': 3,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_preb_4.fastq.gz',
        'tissue': 'bone marrow pre b',
        'biological repetition': 4,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_conv_preb_5.fastq.gz',
        'tissue': 'bone marrow pre b',
        'biological repetition': 5,
        'technical repetition': 1,
        'exposure': 'specific pathogen free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_fo_2_a.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 2,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_gf_fo_2_b.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 2,
        'technical repetition': 2,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_gf_fo_2_c.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 2,
        'technical repetition': 3,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_fo_3_a.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 3,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_fo_3_b.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 3,
        'technical repetition': 2,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_gf_fo_3_c.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 3,
        'technical repetition': 3,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_gf_fo_3_d.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 3,
        'technical repetition': 4,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_fo_4_a.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 4,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_fo_4_b.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 4,
        'technical repetition': 2,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_gf_fo_4_c.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 4,
        'technical repetition': 3,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_gf_fo_4_d.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 4,
        'technical repetition': 4,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_fo_5_a.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 5,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_gf_fo_5_b.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 5,
        'technical repetition': 2,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_gf_fo_5_c.fastq.gz',
        'tissue': 'follicular spleen',
        'biological repetition': 5,
        'technical repetition': 3,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_gf_b1a_2.fastq.gz',
        'tissue': 'peritoneal cavity b1a',
        'biological repetition': 2,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_b1a_3.fastq.gz',
        'tissue': 'peritoneal cavity b1a',
        'biological repetition': 3,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_b1a_4.fastq.gz',
        'tissue': 'peritoneal cavity b1a',
        'biological repetition': 4,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_b1a_5.fastq.gz',
        'tissue': 'peritoneal cavity b1a',
        'biological repetition': 5,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_mz_2.fastq.gz',
        'tissue': 'marginal zone spleen',
        'biological repetition': 2,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_mz_3.fastq.gz',
        'tissue': 'marginal zone spleen',
        'biological repetition': 3,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_mz_4.fastq.gz',
        'tissue': 'marginal zone spleen',
        'biological repetition': 4,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_mz_5.fastq.gz',
        'tissue': 'marginal zone spleen',
        'biological repetition': 5,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_preb_2.fastq.gz',
        'tissue': 'bone marrow pre b',
        'biological repetition': 2,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4G6U/A4G6U l01n01 b6_gf_preb_3.fastq.gz',
        'tissue': 'bone marrow pre b',
        'biological repetition': 3,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_preb_4.fastq.gz',
        'tissue': 'bone marrow pre b',
        'biological repetition': 4,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
    {
        'path': '/volume/albireo/canal/thesis/library/000000000-A4N54/A4N54 l01n01 b6_gf_preb_5.fastq.gz',
        'tissue': 'bone marrow pre b',
        'biological repetition': 5,
        'technical repetition': 1,
        'exposure': 'germ free',
    },
]

def name_library(library):
    buffer = []
    buffer.append('c57bl6')
    buffer.append('b{:0>2}t{:0>2}'.format(library['biological repetition'], library['technical repetition']))
    if library['exposure'] == 'germ free':
        buffer.append('gf')
    else:
        buffer.append('spf')
        
    if library['tissue'] == 'bone marrow pre b':
        buffer.append('preb')

    elif library['tissue'] == 'marginal zone spleen':
        buffer.append('mz')

    elif library['tissue'] == 'peritoneal cavity b1a':
        buffer.append('b1a')
        
    elif library['tissue'] == 'follicular spleen':
        buffer.append('fo')
        
    return ''.join(buffer)

def populate(library):
    name = name_library(library)
    buffer = [
        'gzcat',
        library['path'],
        '|',
        'somatic',
        'populate',
        '--drop',
        '--strain',
        'C57BL/6',
        '--library',
        name,
        '2>',
        '{}/{}.log'.format(LOG, name)
    ]
    for index, literal in enumerate(buffer):
        if ' ' in literal:
            buffer[index] = '"{}"'.format(literal)
    return ' '.join(buffer)

def execute(index):
    if index is not None:
        print(populate(libraries[index]))
    else:
        for library in libraries:
            print(populate(library))
def main():
    index = None
    if len(sys.argv) > 1:
        index = int(sys.argv[1])
    execute(index)

if __name__ == '__main__':
    main()

