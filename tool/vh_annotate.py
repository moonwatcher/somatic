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
import io
import json
from io import StringIO

def to_json(node):
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4)

vh = {}
with io.open('/Users/lg/code/somatic/bootstrap/c57bl6/mouse_ighv_c57bl6.json', 'r') as f:
    document = json.loads(f.read())
    for gene in document:
        vh[gene['id']] = gene
 
with io.open('/Users/lg/Desktop/vh_align.txt', 'r') as file:
    try:
        length = {
            'id': 0,
            'preamble': 0,
            'fr1': 0,
            'cdr1': 0,
            'fr2': 0,
            'cdr2': 0,
            'fr3': 0,
            'cdr3': 0
        }
        for line in file:
            line = line.strip()
            name = line[0:28].strip().strip('>')
            pre = line[28:196].replace('-', '')
            fr1 = line[196:286].replace('-', '')
            cdr1 = line[286:354].replace('-', '')
            fr2 = line[354:419].replace('-', '')
            cdr2 = line[419:463].replace('-', '')
            fr3 = line[463:601].replace('-', '')
            cdr3 = line[601:].replace('-', '')

            if name in vh:
                offset = 0
                node = vh[name]
                node['preamble'] = { 'sequence': { 'nucleotide': pre, 'strand': True } }
                node['fr1'] = { 'sequence': { 'nucleotide': fr1, 'strand': True } }
                node['cdr1'] = { 'sequence': { 'nucleotide': cdr1, 'strand': True } }
                node['fr2'] = { 'sequence': { 'nucleotide': fr2, 'strand': True } }
                node['cdr2'] = { 'sequence': { 'nucleotide': cdr2, 'strand': True } }
                node['fr3'] = { 'sequence': { 'nucleotide': fr3, 'strand': True } }
                node['cdr3'] = { 'sequence': { 'nucleotide': cdr3, 'strand': True } }

                node['preamble']['start'] = offset
                offset += len(pre)
                node['preamble']['end'] = offset
                node['preamble']['accession start'] = node['preamble']['start'] + node['accession start']
                node['preamble']['accession end'] = node['preamble']['end'] + node['accession start']

                node['fr1']['start'] = offset
                offset += len(fr1)
                node['fr1']['end'] = offset
                node['fr1']['accession start'] = node['fr1']['start'] + node['accession start']
                node['fr1']['accession end'] = node['fr1']['end'] + node['accession start']

                node['cdr1']['start'] = offset
                offset += len(cdr1)
                node['cdr1']['end'] = offset
                node['cdr1']['accession start'] = node['cdr1']['start'] + node['accession start']
                node['cdr1']['accession end'] = node['cdr1']['end'] + node['accession start']

                node['fr2']['start'] = offset
                offset += len(fr2)
                node['fr2']['end'] = offset
                node['fr2']['accession start'] = node['fr2']['start'] + node['accession start']
                node['fr2']['accession end'] = node['fr2']['end'] + node['accession start']

                node['cdr2']['start'] = offset
                offset += len(cdr2)
                node['cdr2']['end'] = offset
                node['cdr2']['accession start'] = node['cdr2']['start'] + node['accession start']
                node['cdr2']['accession end'] = node['cdr2']['end'] + node['accession start']

                node['fr3']['start'] = offset
                offset += len(fr3)
                node['fr3']['end'] = offset
                node['fr3']['accession start'] = node['fr3']['start'] + node['accession start']
                node['fr3']['accession end'] = node['fr3']['end'] + node['accession start']

                node['cdr3']['start'] = offset
                offset += len(cdr3)

                # set length
                node['preamble']['length'] = node['preamble']['end'] - node['preamble']['start']
                node['fr1']['length'] = node['fr1']['end'] - node['fr1']['start']
                node['cdr1']['length'] = node['cdr1']['end'] - node['cdr1']['start']
                node['fr2']['length'] = node['fr2']['end'] - node['fr2']['start']
                node['cdr2']['length'] = node['cdr2']['end'] - node['cdr2']['start']
                node['fr3']['length'] = node['fr3']['end'] - node['fr3']['start']

                # update length max
                length['id'] = max(length['id'], len(node['id']))
                length['preamble'] = max(length['preamble'], node['preamble']['length'])
                length['fr1'] = max(length['fr1'], node['fr1']['length'])
                length['cdr1'] = max(length['cdr1'], node['cdr1']['length'])
                length['fr2'] = max(length['fr2'], node['fr2']['length'])
                length['cdr2'] = max(length['cdr2'], node['cdr2']['length'])
                length['fr3'] = max(length['fr3'], node['fr3']['length'])

                # for region in ['fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'preamble']:
                #     if not node[region]['sequence']:
                #         del node[region]

        pattern = '{{: <{}}} {{: <{}}} {{: <{}}} {{: <{}}} {{: <{}}} {{: <{}}} {{: <{}}} {{}}'.format(
            length['id'],
            length['preamble'],
            length['fr1'],
            length['cdr1'],
            length['fr2'],
            length['cdr2'],
            length['fr3']
        )
        # for node in buffer:
        #     print(pattern.format(node['id'], node['preamble']['sequence'], node['fr1']['sequence'], node['cdr1']['sequence'], node['fr2']['sequence'], node['cdr2']['sequence'], node['fr3']['sequence'], node['cdr3']['sequence']))
    except OSError as e:
        print('{} {}'.format(e.strerror, path))

    buffer = vh.values()
    buffer = sorted(buffer, key=lambda x: '' if 'allele' not in x else x['allele'])
    buffer = sorted(buffer, key=lambda x: '' if 'gene' not in x else x['gene'])
    buffer = sorted(buffer, key=lambda x: '' if 'family' not in x else x['family'])
    buffer = sorted(buffer, key=lambda x: '' if 'strain' not in x else x['strain'])

    print(to_json(buffer))


