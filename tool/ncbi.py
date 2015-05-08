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

log = logging.getLogger('Accession')

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
                            node = INSDSeq
                            break
    return node

def to_json(node):
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4)

print(to_json(fetch_from_ncbi(sys.argv[1])))






