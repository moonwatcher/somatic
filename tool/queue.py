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
import io

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


libraries = []
with io.open('/Users/lg/code/somatic/bootstrap/library.json', 'r') as f:
    document = json.loads(f.read())
    document = sorted(document, key=lambda x: 0 if 'id' not in x['head'] else x['head']['id'])
    document = sorted(document, key=lambda x: 0 if 'technical repetition' not in x['body'] else x['body']['technical repetition'])
    document = sorted(document, key=lambda x: 0 if 'biological repetition' not in x['body'] else x['body']['biological repetition'])
    document = sorted(document, key=lambda x: '' if 'tissue' not in x['body'] else x['body']['tissue'])
    document = sorted(document, key=lambda x: '' if 'exposure' not in x['body'] else x['body']['exposure'])

    for library in document:
        if 'reference' not in library['body']:
            if library['head']['tissue'] in [
                'immature adult bone marrow',
                'immature adult spleen',
                'immature fetal liver',
                'pre b adult bone marrow',
                'pre b fetal liver'
            ]:
                libraries.append(library)


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
    buffer = [
        'gzcat',
        library['body']['path'],
        '|',
        'somatic',
        'populate',
        '--drop',
        '--strain',
        'C57BL/6',
        '--library',
        library['head']['id'],
        '2>',
        '{}/{}.log'.format(LOG, library['head']['id'])
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
    # record()

