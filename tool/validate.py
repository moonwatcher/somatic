#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import logging
import re
import json
import uuid
import io
import hashlib
import math

from io import StringIO, BytesIO
from datetime import timedelta, datetime
from argparse import ArgumentParser
from subprocess import Popen, PIPE

from bson.objectid import ObjectId
from pymongo.son_manipulator import SONManipulator
from pymongo import MongoClient, DESCENDING, ASCENDING
from pymongo.errors import BulkWriteError

import xmltodict
import urllib.request, urllib.parse, urllib.error
import urllib.parse
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from http.client import BadStatusLine

log = logging.getLogger('Validate')

def reverse(sequence):
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
    return ''.join([complement[b] for b in list(sequence.strip())][::-1])

def fetch_from_ncbi(id):
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

    def fetch(url):
        content = None
        request = Request(url, None, { 'Accept': 'application/xml' })
        
        try:
            response = urlopen(request)
        except BadStatusLine as e:
            log.warning('Bad http status error when requesting %s', url)
        except HTTPError as e:
            log.warning('Server returned an error when requesting %s: %s', url, e.code)
        except URLError as e:
            log.warning('Could not reach server when requesting %s: %s', url, e.reason)
        else:
            content = StringIO(response.read().decode('utf8'))
            if content.read(22) == 'Nothing has been found':
                content = None
            else:
                content.seek(0)
        return parse(content)

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
            document = xmltodict.parse(content.getvalue())
            document = normalize_document(document, transform)
        return document

    node = None
    if id is not None:
        url = 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&dopt=gbc_xml&val={}'
        document = fetch(url.format(id))
        if 'INSDSet' in document:
            for INSDSet in document['INSDSet']:
                if 'INSDSeq' in INSDSet:
                    for INSDSeq in INSDSet['INSDSeq']:
                        if 'INSDSeq_moltype' in INSDSeq and INSDSeq['INSDSeq_moltype'] in ['DNA', 'mRNA', 'RNA']:
                            if 'INSDSeq_sequence' in INSDSeq:
                                INSDSeq['INSDSeq_sequence'] = INSDSeq['INSDSeq_sequence'].upper()
                            node = INSDSeq
                            break
    return node

def to_json(node):
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4)

def search(space, query):
    result = []
    start = space.find(query)
    while start >= 0 and start < len(space):
        end = start + len(query)
        result.append([start, end])
        start = space.find(query, start + 1)
    return result

def validate(gene):
    if 'accession' in gene and 'start' in gene and 'end' in gene and 'strand' in gene:
        if gene['accession'] in documents:
            document = documents[gene['accession']]
        else:
            document = fetch_from_ncbi(gene['accession'])
            documents[gene['accession']] = document

        r = document['INSDSeq_sequence'][gene['start']:gene['end']].upper()
        if not gene['strand']:
            r = reverse(r)

        if r != gene['sequence']:
            complement(gene, True)
            r = document['INSDSeq_sequence'][gene['start']:gene['end']].upper()
            if r != gene['sequence']:
                print('ERROR: {} is invalid'.format(gene['id']))
        else:
            print('{} is valid'.format(gene['id']))
    else:
        print('ERROR: {} is missing metadata'.format(gene['id']))

def complement(gene, force=False):
    if force or ('start' not in gene or 'end' not in gene):
        if gene['accession'] in documents:
            document = documents[gene['accession']]
        else:
            document = fetch_from_ncbi(gene['accession'])
            documents[gene['accession']] = document
        result = search(document['INSDSeq_sequence'], gene['sequence'])
        if len(result) == 1:
          r = result[0]
          print('new position found for {} in {} : {}'.format(gene['id'], gene['accession'], r))
          gene['start'] = r[0]
          gene['end'] = r[1]
        elif len(result) == 0:
          print('{} not found in {}'.format(gene['id'], gene['accession']))
          print(document['INSDSeq_sequence'])
          print(gene['sequence'])
        else:
          print('multiple positions found for {} in {} : {}'.format(gene['id'], gene['accession'], result))

BN000872 = set([
    "J558.89pg.195",
    "J558.88.194",
    "J558.87.193",
    "J558.86.192",
    "J558.85.191",
    "J558.84.190",
    "J558.83.189",
    "J558.82pg.188",
    "J558.81.187",
    "J558.80.186",
    "3609.16pg.185",
    "J558.79.184",
    "3609.15pg.183",
    "J558.78.182",
    "3609.14pg.181",
    "J558.77.180",
    "J558.76pg.179",
    "3609.13pg.178",
    "J558.75.177",
    "J558.74.176",
    "J558.73pg.175",
    "3609.12.174",
    "J558.72.173",
    "J558.71pg.172",
    "J558.70pg.171",
    "J558.69.170",
    "3609.11.169",
    "J558.68pg.168",
    "3609.10pg.167",
    "J558.67.166",
    "J558.66.165",
    "3609.9.164",
    "J558.65.163",
    "J558.64.162",
    "PG.19.161",
    "3609.8pg.160",
    "J558.63pg.159",
    "J558.62pg.158",
    "J558.61.157",
    "J558.60pg.156",
    "J558.59.155",
    "J558.58.154",
    "3609.7.153",
    "J558.57pg.152",
    "3609.6pg.151",
    "J558.56.150",
    "J558.55.149",
    "J558.54.148",
    "3609.5.147",
    "J558.53.146",
    "J558.52.145",
    "J558.51pg.144",
    "J558.50.143",
    "3609.4.142",
    "J558.49.141",
    "J558.48pg.140",
    "3609.3.139",
    "3609.2pg.138",
    "J558.47.137",
    "J558.46pg.136",
    "J558.45pg.135",
    "J558.44pg.134",
    "J558.43.133",
    "J558.42.132",
    "J558.41pg.131",
    "J558.40pg.130",
    "J558.39.129",
    "J558.38pg.128",
    "J558.37.127",
    "J558.36.126",
    "J558.35pg.125",
    "J558.34.124",
    "J558.33pg.123",
    "J558.32pg.122",
    "J558.31.121",
    "J558.30pg.120",
    "J558.29pg.119",
    "J558.28pg.118",
    "J558.27pg.117",
    "J558.26.116",
    "J558.25pg.115",
    "J558.24pg.114",
    "J558.23.113",
    "J558.22.112",
    "J558.21pg.111",
    "J558.20pg.110",
    "J558.19.109",
    "J558.18.108",
    "J558.17pg.107",
    "J558.16.106",
    "J558.15pg.105",
    "J558.14pg.104",
    "J558.13.103",
    "J558.12.102",
    "J558.11pg.101",
    "J558.10pg.100",
    "J558.9.99",
    "J558.8.98",
    "J558.7pg.97",
    "J558.6.96",
    "VH15.1.95",
    "J558.5pg.94",
    "J558.4.93",
    "PG.18.92",
    "VH10.3.91",
    "J558.3.90",
    "VH10.2pg.89",
    "J558.2.88",
    "PG.17.87",
    "VH10.1.86",
    "J558.1.85",
    "3609.1.84",
    "J606.5.83",
    "J606.4.82",
    "J606.3.81",
    "J606.2.80",
    "J606.1.79",
    "VH12.1.78",
    "3609N.2.77",
    "PG.16.76",
    "PG.15.75",
    "36-60.8.74",
    "PG.14.73",
    "36-60.7pg.72",
    "VGAM3.8-4-71",
    "36-60.6.70",
    "PG.13.69",
    "3609N.1pg.68",
    "36-60.5.67",
    "36-60.4.66",
    "S107.4.65",
    "36-60.3.64",
    "SM7.4.63",
    "S107.3.62",
    "VGAM3.8-3-61",
    "PG.12.60",
    "VGAM3.8-2-59",
    "PG.11.58",
    "VGAM3.8-1-57",
    "PG.10.56",
    "VH16.1.55",
    "SM7.3.54",
    "VH11.2.53",
    "PG.9.52",
    "36-60.2pg.51",
    "X24.2.50",
    "SM7.2.49",
    "VH11.1.48",
    "PG.8.47",
    "36-60.1.46",
    "X24.1pg.45",
    "SM7.1.44",
    "S107.2pg.43",
    "S107.1.42",
    "PG.7.41",
    "Q52.13.40",
    "Q52.12pg.39",
    "7183.21pg.38",
    "7183.20.37",
    "7183.19.36",
    "7183.18.35",
    "Q52.11.34",
    "Q52.10.33",
    "PG.6.32",
    "7183.17pg.31",
    "PG.5.30",
    "Q52.9.29",
    "PG.4.28",
    "7183.16.27",
    "7183.15pg.26",
    "7183.14.25",
    "7183.13pg.24",
    "PG.3.23",
    "Q52.8.22",
    "PG.2.21",
    "7183.12.20",
    "7183.11pg.19",
    "Q52.7.18",
    "Q52.6pg.17",
    "7183.10pg.16",
    "7183.9.15",
    "7183.8pg.14",
    "Q52.5.13",
    "Q52.4pg.12",
    "PG.1.11",
    "7183.7.10",
    "7183.6pg.9",
    "Q52.3.8",
    "7183.5pg.7",
    "7183.4.6",
    "7183.3pg.5",
    "Q52.2.4",
    "7183.2.3",
    "Q52.1pg.2",
    "7183.1pg.1",
])
lookup = {}

documents = {}
# for k,region in genes.items():
#     for gene in region:
#         validate(gene)
#document = json.loads(sys.stdin.read())

notconfirmed = set([
    'IGHV3-6*03', 
    'IGHV3S7*01', 
    'IGHV10S4*01', 
    'IGHV10S5*01', 
    'IGHV12S2*01', 
    'IGHV13-2*03', 
    'IGHV14S4*01', 
    'IGHV15-2*02'
])


notidentified = set([
    'IGHV1-18*02',
    'IGHV1-18*03',
    'IGHV1-33*02',
    'IGHV1-42*02',
    'IGHV1-42*03',
    'IGHV1-53*02',
    'IGHV1-53*03',
    'IGHV1-53*04',
    'IGHV1-54*03',
    'IGHV1-55*02',
    'IGHV1-55*03',
    'IGHV1-55*04',
    'IGHV1-56*02',
    'IGHV1-62-3*02',
    'IGHV1-64*02',
    'IGHV1-69*03',
    'IGHV1-72*02',
    'IGHV1-72*03',
    'IGHV1-72*05',
    'IGHV1-74*02',
    'IGHV1-74*03',
    'IGHV1-84*02',
    'IGHV1S100*01',
    'IGHV1S101*01',
    'IGHV1S101*02',
    'IGHV1S103*01',
    'IGHV1S107*01',
    'IGHV1S108*01',
    'IGHV1S110*01',
    'IGHV1S111*01',
    'IGHV1S112*01',
    'IGHV1S112*02',
    'IGHV1S113*01',
    'IGHV1S113*02',
    'IGHV1S118*01',
    'IGHV1S120*01',
    'IGHV1S120*02',
    'IGHV1S121*01',
    'IGHV1S122*01',
    'IGHV1S124*01',
    'IGHV1S126*01',
    'IGHV1S127*01',
    'IGHV1S130*01',
    'IGHV1S132*01',
    'IGHV1S134*01',
    'IGHV1S135*01',
    'IGHV1S136*01',
    'IGHV1S137*01',
    'IGHV1S20*02',
    'IGHV1S21*02',
    'IGHV1S36*02',
    'IGHV1S44*01',
    'IGHV1S46*01',
    'IGHV1S65*01',
    'IGHV1S65*02',
    'IGHV1S65*03',
    'IGHV1S67*01',
    'IGHV1S67*02',
    'IGHV1S68*01',
    'IGHV1S68*02',
    'IGHV1S70*01',
    'IGHV1S72*01',
    'IGHV1S73*01',
    'IGHV1S74*01',
    'IGHV1S75*01',
    'IGHV1S75*02',
    'IGHV1S78*01',
    'IGHV1S81*01',
    'IGHV1S81*02',
    'IGHV1S82*01',
    'IGHV1S83*01',
    'IGHV1S84*01',
    'IGHV1S87*01',
    'IGHV1S88*01',
    'IGHV1S92*01',
    'IGHV1S95*01',
    'IGHV1S96*01',
    'IGHV5-12*03',
    'IGHV5-15*03',
    'IGHV5-15*04',
    'IGHV5-16*02',
    'IGHV5-17*03',
    'IGHV5-2*03',
    'IGHV5-4*03',
    'IGHV5-6*02',
    'IGHV5-6*03',
    'IGHV5-9*04',
    'IGHV5S12*01',
    'IGHV5S21*01',
    'IGHV5S24*01',
    'IGHV6S2*01',
    'IGHV6S3*01',
    'IGHV6S4*01',
    'IGHV8-5*02',
    'IGHV8-8*02',
    'IGHV8-9*02',
    'IGHV8-9*03',
    'IGHV8S6*01',
    'IGHV8S9*01',
])

original = None
with io.open('/Users/lg/code/somatic/bootstrap/mouse_ighv.json', 'r') as f:
    original = json.loads(f.read())

for gene in original:
    lookup[gene['id']] = gene

buffer = []
for gene in lookup.values():
    gene['confirmed'] = True
    gene['identified'] = True
    gene['verified'] = True

    # if 'name' in gene: del gene['name']

    if 'imgt allele name' in gene:
        if gene['imgt allele name'] in notconfirmed:
            print(gene['id'])
            gene['confirmed'] = False

        if gene['imgt allele name'] in notidentified:
            print(gene['id'])
            gene['identified'] = False

        if not gene['identified'] or not gene['confirmed']:
            print(gene['id'])
            gene['verified'] = False

# document = json.loads(sys.stdin.read())

# family = {
#     'imgt': {
#         'IGHV1':'J558',
#         'IGHV2':'Q52',
#         'IGHV3':'36-60',
#         'IGHV4':'X24',
#         'IGHV5':'7183',
#         'IGHV6':'J606',
#         'IGHV7':'S107',
#         'IGHV8':'3609v',
#         'IGHV9':'VGAM3.8',
#         'IGHV10':'VH10',
#         'IGHV11':'VH11',
#         'IGHV12':'VH12',
#         'IGHV13':'3609N',
#         'IGHV14':'IGHV14',
#         'IGHV15':'IGHV15',
#         'IGHV16':'IGHV16',
#         'IGHV(II)':'VHPS2',
#         'IGHV(III)':'VHPS3',
#     }
# }

# family['lit'] = {}
# for k,v in family['imgt'].items():
#     family['lit'][v] = k

# buffer = []
# for gene in document:
#     if 'family name' in gene:
#         gene['family'] = gene['family name']
#         del gene['family name']


#     if 'imgt family name' in gene and 'family' not in gene:
#         gene['family'] = family['imgt'][gene['imgt family name']]
#     #     print('adding family name {} to {}'.format(gene['family name'], gene['id']))

#     if 'family' in gene and 'imgt family name' not in gene:
#         gene['imgt family name'] = family['lit'][gene['family']]
#     #     print('adding imgt family name {} to {}'.format(gene['imgt family name'], gene['id']))

#     x = [ 
#         gene['imgt family name'],
#         gene['family'],
#         gene['name'],
#         gene['id'],
#     ]
#     x.append('None' if 'strain' not in gene else gene['strain'])
#     x.append('None' if 'imgt allele name' not in gene else gene['imgt allele name'])
#     x.append(gene['accession'])
#     x.append(gene['start'])
#     x.append(gene['end'])
#     buffer.append(x)

buffer = sorted(buffer, key=lambda x: x['gene'])
buffer = sorted(buffer, key=lambda x: x['family'])

#print('\n'.join([ ','.join([str(i) for i in x]) for x in buffer ]))

print(to_json(buffer))




