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

document = None
with io.open('../bootstrap/c57bl6/library_c57bl6.json', 'r') as f:
    document = json.loads(f.read())
    document = sorted(document, key=lambda x: 0 if 'id' not in x['head'] else x['head']['id'])
    document = sorted(document, key=lambda x: 0 if 'technical repetition' not in x['body'] else x['body']['technical repetition'])
    document = sorted(document, key=lambda x: 0 if 'biological repetition' not in x['body'] else x['body']['biological repetition'])
    document = sorted(document, key=lambda x: '' if 'tissue' not in x['body'] else x['body']['tissue'])
    document = sorted(document, key=lambda x: '' if 'exposure' not in x['body'] else x['body']['exposure'])

print('{} | {} | {} | {} | {} | {}'.format('name', 'exposure', 'tissue', 'B', 'T', 'Reference'))
print('{} | {} | {} | {} | {} | {}'.format('-' * (len('name') + 2), '-' * (len('exposure') + 2), '-' * (len('tissue') + 2), '-' * (len('B') + 2), '-' * (len('T') + 2), '-' * (len('Reference') + 2)))
for library in document:
    if 'reference' not in library['body']:
        print('{} | {} | {} | {} | {} |'.format(
            library['head']['id'],
            '' if 'exposure' not in library['body'] else library['body']['exposure'],
            '' if 'tissue' not in library['body'] else library['body']['tissue'],
            '' if 'biological repetition' not in library['body'] else library['body']['biological repetition'],
            '' if 'technical repetition' not in library['body'] else library['body']['technical repetition']
            )
        )
    else:
        print('{} | {} | {} | | | {}'.format(
            library['head']['id'],
            '' if 'exposure' not in library['head'] else library['head']['exposure'],
            '' if 'tissue' not in library['head'] else library['head']['tissue'],
            ', '.join(library['body']['reference'])
            )
        )



