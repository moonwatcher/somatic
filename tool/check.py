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

def to_json(node):
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4)

index = {
    'good': {},
    'bad': {}
}

good = [
    {
        "accession": "BN000872",
        "allele": "36-60.1.46*01",
        "confirmed": True,
        "end": 2106764,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.1.46",
        "id": "36-60.1.46*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-1*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-1",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GATGTGCAGCTTCAGGAGTCAGGACCTGGCATGGTGAAACCTTCTCAGTCACTTTCCCTCACCTGCACTGTCACTGGCTACTCCATCACCAGTGGTTATGACTGGCACTGGATCCGACATTTTCCAGGAAACAAACTGGAGTGGATGGGCTACATAAGCTACAGTGGTAGCACTAACTACAACCCATCCCTCAAAAGTCGAATCTCCATCACTCATGACACATCTAAGAACCATTTCTTCCTGAAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2106468,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "36-60.2pg.51*01",
        "confirmed": True,
        "end": 2037309,
        "family": "36-60",
        "framed": True,
        "functionality": "P",
        "gene": "36-60.2pg.51",
        "id": "36-60.2pg.51*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-2*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-2",
        "length": 295,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GATGTGCAGCTTCAGGAGTCAGGACCTGGCCTGGTGAAACCTTCTCAGTCACTTTCCCTCACCTGCACTGTCACTGGCTACTCCATCACCAGTGGTTATGACTGGAACTGGATCCGGCAGTTTCCAGGAAACAAACTGGAGTGGATGGGCTACATAAGCTACAGTGTAGCAAGTACTACAACCCATCTCTCAAAAGTCGAATCTCTATCACTCGAGACACATCCAAGAACCAGTTCTCCCTGGAATTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2037014,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "36-60.3.64*01",
        "confirmed": True,
        "end": 1874673,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.3.64",
        "id": "36-60.3.64*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-3*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-3",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GATGTGCAGCTTCAGGAGTCAGGACCTAGCCTGGTGAGACCTTCTCAGACACTCTCCCTTACCTGCACTGTTACTGGCTTCTCCATCAACAGTGATTGTTACTGGATCTGGATCCGGCAGTTTCCAGGAAACAAACTGGAGTACATCGGGTACACATTCTACAGTGGTATCACTTACTACAACCCATCTCTTGAAAGTCGAACGTACATAACGCGTGACACATCTAAGAACCAGTTCTCACTGAAGTTGAGTTCTGTGACTACTGAGGACACAGCCACTTACTACTGTGCGAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1874377,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "36-60.4.66*01",
        "confirmed": True,
        "end": 1817495,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.4.66",
        "id": "36-60.4.66*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-4*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-4",
        "length": 299,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GATGTACAGCTTCAGGAGTCAGGACCTGCCCTGGTGAAGCCTTCTCAGACAGTGTCCCTCACCTGCACTGTCACTGGCTACTCTATCACTAATGGTAATCACTGGTGGAACTGGATCCGGCAGGTTTCAGGAAGCAAACTGGAGTGGATAGGGTACATAAGCTCCAGTGGTAGCACTGACAGCAATCCATCTCTCAAAAGTCGAATCTCCATCACTAGAGACACTTCCAAGAACCAGTTATTCCTGCAGTTGAACTCTGTGACTACTGAAGATATAGCCACATATTACTGTGCAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1817196,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "36-60.5.67*01",
        "confirmed": True,
        "end": 1808460,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.5.67",
        "id": "36-60.5.67*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-5*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-5",
        "length": 299,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GATGTGCAGCTTCAGGAGTCAGGACCTGGCCTGGTGAAACCTTCTCAGACAGTGTTCCTCACCTGCACTGTCACTGGCATTTCCATCACCACTGGAAATTACAGGTGGAGCTGGATCCGGCAGTTTCCAGGAAACAAACTGGAGTGGATAGGGTACATATACTACAGTGGTACCATTACCTACAATCCATCTCTCACAAGTCGAACCACCATCACTAGAGACACTCCCAAGAACCAGTTCTTCCTGGAAATGAACTCTTTGACTGCTGAGGACACAGCCACATACTACTGTGCACGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1808161,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "36-60.6.70*01",
        "confirmed": True,
        "end": 1782958,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.6.70",
        "id": "36-60.6.70*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-6*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-6",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GATGTACAGCTTCAGGAGTCAGGACCTGGCCTCGTGAAACCTTCTCAGTCTCTGTCTCTCACCTGCTCTGTCACTGGCTACTCCATCACCAGTGGTTATTACTGGAACTGGATCCGGCAGTTTCCAGGAAACAAACTGGAATGGATGGGCTACATAAGCTACGATGGTAGCAATAACTACAACCCATCTCTCAAAAATCGAATCTCCATCACTCGTGACACATCTAAGAACCAGTTTTTCCTGAAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1782662,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "36-60.7pg.72*01",
        "confirmed": True,
        "end": 1756925,
        "family": "36-60",
        "framed": True,
        "functionality": "P",
        "gene": "36-60.7pg.72",
        "id": "36-60.7pg.72*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-7*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-7",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTACAGCTCCAGGAGTCGAGACCTAGCCTGGTGAAACCTTCTCAGTCATTGTCCTTCACCTGCTCTGTCACTGGCTACTCCATCACCAGTGGTTATGACTGGAACTGGATCCGGCAGTTTCCAGGAAACAAACTGTAATGGATGAGGTACATACACAAGAGTGGTAGCACTAACTACAATCCATCTCTCAAAAGTCAAGTCTCCATCACTCGAGACACTTCCAAGAACCAGTTCTACCTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCATATATTACTGTGCAAACGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1756629,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "36-60.8.74*01",
        "confirmed": True,
        "end": 1748735,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.8.74",
        "id": "36-60.8.74*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-8*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-8",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGCAGCTTCAGGAGTCAGGACCTGGCCTGGCAAAACCTTCTCAGACTCTGTCCCTCACCTGTTCTGTCACTGGCTACTCCATCACCAGTGATTACTGGAACTGGATCCGGAAATTCCCAGGGAATAAACTTGAGTACATGGGGTACATAAGCTACAGTGGTAGCACTTACTACAATCCATCTCTCAAAAGTCGAATCTCCATAACTCGAGACACATCCAAGAACCAGTATTACCTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGATA",
            "read frame": 0,
            "strand": True
        },
        "start": 1748442,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609N.1pg.68*01",
        "confirmed": True,
        "end": 1803509,
        "family": "3609N",
        "framed": True,
        "functionality": "P",
        "gene": "3609N.1pg.68",
        "id": "3609N.1pg.68*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV13-1*01",
        "imgt family name": "IGHV13",
        "imgt gene name": "IGHV13-1",
        "length": 304,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTACAGCTGGTAGAGACAGGAGGAGGCTTGGTGCAGCCTGGAAACTCTCTAAAACTTTCCTGTGCCACTTCGGGATACCCTTTTTATGACTACTGGATGGATTGGGTCCGCCATTCTCCAGAAAAGGGGCTGGAGTGGGTTGCTCGAATTGCAACAAAAACTCATAATTATGCAACGTACTATGCAGAGTCTGTGAAAGGCCGATTCATCGTCTCAAGAGATGATTCCAAAAGCAGTGCATACATGCAGATGAACAGCTTAAGAAAGGAAGACACTGCCGTTTATTACTGTGCAAGAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1803205,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609N.2.77*01",
        "confirmed": True,
        "end": 1713348,
        "family": "3609N",
        "framed": True,
        "functionality": "F",
        "gene": "3609N.2.77",
        "id": "3609N.2.77*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV13-2*01",
        "imgt family name": "IGHV13",
        "imgt gene name": "IGHV13-2",
        "length": 300,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTGCAGCTTGTAGAGACCGGGGGAGGCTTGGTGAGGCCTGGAAATTCTCTGAAACTCTCCTGTGTTACCTCGGGATTCACTTTCAGTAACTACCGGATGCACTGGCTTCGCCAGCCTCCAGGGAAGAGGCTGGAGTGGATTGCTGTAATTACAGTCAAATCTGATAATTATGGAGCAAATTATGCAGAGTCTGTGAAAGGCAGATTCGCCATTTCAAGAGATGATTCAAAAAGCAGTGTCTACCTAGAGATGAACAGATTAAGAGAGGAAGACACTGCCACTTATTTTTGTAGTAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1713048,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.1.84*01",
        "confirmed": True,
        "end": 1608819,
        "family": "3609v",
        "framed": True,
        "functionality": "O",
        "gene": "3609.1.84",
        "id": "3609.1.84*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-2*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-2",
        "length": 300,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTACTCTGAAAGTGTCTGGCCCTGGGATATTGCAGCCATCACAGACTCTCAGCCTGGCCTGTACTTTCTCTGGGATTTCACTGAGTACTTCTGGTATGGGTTTGAGCTGGCTTCGTAAGCCCTCAGGGAAGGCTTTAGAGTGGCTGGCAAGCATTTGGAATAATGATAACTACTACAACCCATCTTTGAAGAGCCGGCTCACAATCTCCAAGGAGACCTCCAACTACCAAGTATTCCTTAAACTCACCAGTGTGGACACTGCAGATTCTGCCACATACTACGGTGCTTGGAGAGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 1608519,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.10pg.167*01",
        "confirmed": True,
        "end": 551738,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.10pg.167",
        "id": "3609.10pg.167*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-10*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-10",
        "length": 301,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAAGTTACGCTTTAACAGTCTGGCCCTGGGAAATTGCAGCCCTCCCCGACCTTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTATGGAATGATGGTGAGCTGGATGTGTCAGCCTTCAGGGAAGGGTGTGGAGTGGCTGGCACACATTTGGTGGAATGATGATAAGGGCTACAACCCATCCCTGAGGAGCCGGCTGACAATCTCCAAAGATACCCCCAACACCCAGGTATTCCTTAAGATCTCCAGTGTGGTCACTGCAGATACTACCACATACTACTGTGCTTGAGGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 551437,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.11.169*01",
        "confirmed": True,
        "end": 503929,
        "family": "3609v",
        "framed": True,
        "functionality": "F",
        "gene": "3609.11.169",
        "id": "3609.11.169*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-11*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-11",
        "length": 301,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATTACTCAGAAAGAGTCTGGCCCTGGGATATTGCAGCCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTCTGGTATGGGTGTAGGCTGGATTCATCAGCCTTCAGGGAATGGTCTGGAGTGGCTGGCACACATTTGGTGGAATGATAATAAGTACTATAACACAGCCCTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCGCCAGTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAATAG",
            "read frame": 0,
            "strand": True
        },
        "start": 503628,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.12.174*01",
        "confirmed": True,
        "end": 423133,
        "family": "3609v",
        "framed": True,
        "functionality": "F",
        "gene": "3609.12.174",
        "id": "3609.12.174*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-12*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-12",
        "length": 301,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTACTCTGAAAGAGTCTGGCCCTGGGATATTGCAGTCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTCTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGAAAGGGTCTGGAGTGGCTGGCACACATTTACTGGGATGATGACAAGCGCTATAACCCATCCCTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAGAAACCAGGTATTCCTCAAGATCACCAGTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAAGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 422832,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.13pg.178*01",
        "confirmed": True,
        "end": 305743,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.13pg.178",
        "id": "3609.13pg.178*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-13*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-13",
        "length": 301,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTACTCTGAAAGAGTCTGGCCCTGGTATATTGCAGCCCTCCCAGACCCTCAGTCTGACCTGTTCTTTCTCTGGGTTTTCACTGAGCACTTTTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTGGAGTGGCTGGCACACATTTATTGGGATGATGACAAGCACTATAACCCATCCTTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCACTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAAGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 305442,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.14pg.181*01",
        "confirmed": True,
        "end": 263628,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.14pg.181",
        "id": "3609.14pg.181*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-14*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-14",
        "length": 301,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAAGTTACTCTAAAAGAGTATGGCCCTGGGAAACTGTAGCCCTCTCAGACCTTCAGTCTGACTTGTACTTTCTCTGGGTTTTCACTGAGCACTTATGGAATGATGGTGAGCTGCATGTGTCAGCCTTCAGGGAAGGGTCTGGTGTGGCTGGCACTCATTTGGTGGAATAATGATAAGGGCTACAACCCATTCCTGAGGAGCCAGCTGACAATCTCCAAGGATACATCCAACAACCAGGTATTCCTCAAGATCACCAGTGTGGACCCTGCAGATACTGCCACATACTACTGTGCTTGAGGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 263327,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.15pg.183*01",
        "confirmed": True,
        "end": 228044,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.15pg.183",
        "id": "3609.15pg.183*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-15*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-15",
        "length": 299,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGTTTACTCTAAAAGAGTCTGGCCCTGGGATACTGCAGACCTCCAAGACGCTCAGTCTTACTTGCTCTTCCTTGGGTTTTCACAGAGAACTTCTGGTTTGGTCATGAGCGGGATCTGTCAGCCATCAGGGAACTGTCTGAGGTAGTTGGCATTAGTTGATGGAATGATGACAAGTACTACAACCCTTCTCTGAAGAGCTGGATCACAATATCCAAGGACACCTCCACCAACCAACCATTCCTCAAGACCACCACTGTGGACACTGCAGACACTGCCACATATTAATGTGCTCCAGGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 227745,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.16pg.185*01",
        "confirmed": True,
        "end": 214051,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.16pg.185",
        "id": "3609.16pg.185*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-16*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-16",
        "length": 299,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGTTTACTCTAAAAGAGTCTGGCCCTGGGATACTGCAGACCTCCAAGACGCTCAGTCTTACTTGCTCTTCCTTGGGTTTTCACAGAGAACTTCTGGTTTGGTCATGAGTGGGATCTGTCAGCCATCAGGGAACTGTCTGAGGTAGTTGGCATTAATTGATGGAATGATGATAAGTACTACAACCCTTCTCTGAAGAGCTGGATCACAATATCCAAGGACACCTCCACCAACCAACCATTCCTCAAGACCACCACTGTGGACACTGCAGATAGTGCCACATATTAATGTGCTCCAGGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 213752,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.2pg.138*01",
        "confirmed": True,
        "end": 1067659,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.2pg.138",
        "id": "3609.2pg.138*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-3*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-3",
        "length": 295,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTCTGTCTTAGGTTACTTTGAAAGTGTCTGGCCCTGGCGTATTGCATCCCTCCCAGGCTCACAGTCTGACTTGCTCATTGTGTTTTCACTGAGCACTTCTGGTAGGACCATGATTCATCAGCCTTCAGGGAAAGGTTTGACAACAATTTGTTGGAATTATGATAAGTGCTACAACACATCTATTAAGAGTCAGCTCACAGTCTCCAAGGACACCTCAAACAACCAAGCATTACTCCAGATCATCAATGTGGACACTGCAAATACTGACACATACTACGGTGCCTCAGTTTT",
            "read frame": 0,
            "strand": True
        },
        "start": 1067364,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.3.139*01",
        "confirmed": True,
        "end": 1047058,
        "family": "3609v",
        "framed": True,
        "functionality": "O",
        "gene": "3609.3.139",
        "id": "3609.3.139*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-4*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-4",
        "length": 301,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATTACTCTGAAACAGTCTGGCCCTGGGATAGTGCAGCCATCCCAGCCCGTCAGACTTACTTGCACTTTCTCTGGGTTTTCACTGAGCACTTCTGGTATAGGTGTAACCTGGATTCGTCAGCCCTCAGGGAAGGGTCTGGAGTGGCTGGCAACCATTTGGTGGGATGATGATAACCGCTACAACCCATCTCTAAAGAGCAGGCTCGCAGTCTCCAAAGACACCTCCAACAACCAAGCATTCCTGAATATCATCACTGTGGAAACTGCAGATACTGCCATATACTACTGTGCTCAGAGTG",
            "read frame": 0,
            "strand": True
        },
        "start": 1046757,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.4.142*01",
        "confirmed": True,
        "end": 1003518,
        "family": "3609v",
        "framed": True,
        "functionality": "F",
        "gene": "3609.4.142",
        "id": "3609.4.142*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-5*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-5",
        "length": 301,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTACTCTGAAAGAGTCTGGCCCTGGGATATTGCAGCCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTCTAATATGGGTATAGGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTAGAGTGGCTGGCACACATTTGGTGGAATGATGATAAGTACTATAACCCATCCCTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCAGTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCAAATAG",
            "read frame": 0,
            "strand": True
        },
        "start": 1003217,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.5.147*01",
        "confirmed": True,
        "end": 905301,
        "family": "3609v",
        "framed": True,
        "functionality": "F",
        "gene": "3609.5.147",
        "id": "3609.5.147*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-6*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-6",
        "length": 301,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTACTCTGAAAGAGTCTGGCCCTGGTATATTGCAGCCCTCGCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGTACTTTTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGGAAGGATCTGGAGTGGCTGGCACACATTTATTGGGATGATGACAAGCACTATAACCCATCCTTGAAGAGCCAGCTCAGAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCACTGTGGACACTGTAGATACTGCCACATACTACTGTGCTCGAAGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 905000,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.6pg.151*01",
        "confirmed": True,
        "end": 812941,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.6pg.151",
        "id": "3609.6pg.151*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-7*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-7",
        "length": 292,

        "organism name": "mus musculus",
        "read frame": 2,
        "region": "VH",
        "sequence": {
            "nucleotide": "TTATTCTGAAAGAGTCTGGCCCTGGAATATTGCAGCCGTCCCAGACCCTAAGTCTGACTTGTTCTTTCTCTGGGTTTTCATTTAGGACTTAAGGTATGGCAGTGAGCTGGATGCGTCAGCCTTTAGGGAAGGGTTTGGAGTGGCTGGCACAAATTGGGTCAGATGATAGTAAGCTCTATAACCCATCTCTGAAGAGTCGAACCACAATCTCCAAGGATACCTCCAACAACCATGTTTTCCTCAAGATCACCAGTGAGGACACTGAAGATTCTGCCACATACTACTGTGCTCA",
            "read frame": 2,
            "strand": True
        },
        "start": 812649,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.7.153*01",
        "confirmed": True,
        "end": 777016,
        "family": "3609v",
        "framed": True,
        "functionality": "F",
        "gene": "3609.7.153",
        "id": "3609.7.153*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-8*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-8",
        "length": 301,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTACTCTGAAAGAGTCTGGCCCTGGGATATTGCAGCCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTTTGGTATGGGTGTAGGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTGGAGTGGCTGGCACACATTTGGTGGGATGATGATAAGTACTATAACCCAGCCCTGAAGAGTCGGCTCACAATCTCCAAGGATACCTCCAAAAACCAGGTATTCCTCAAGATCGCCAATGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAATAG",
            "read frame": 0,
            "strand": True
        },
        "start": 776715,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.8pg.160*01",
        "confirmed": True,
        "end": 666661,
        "family": "3609v",
        "framed": True,
        "functionality": "O",
        "gene": "3609.8pg.160",
        "id": "3609.8pg.160*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-8-1*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-8-1",
        "length": 306,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTATTCTGAAAGAGTCTGGCCCTGGAATATTGCAGCCCTCCCAGACCCTAAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTTAGCACTTATGGTATGACAGTGAGCTGGATGCGTCAGCCTTCAAGGAAGGGTTTGGAGTGGCTGGCACAAATTGGGTCAGATGATAGTAAGCTCTATAACCCATCTCTGAAGAGCCGAATCACAATCTCCAAGGAAACCTCCAACAACCATGTATTCCTCAAGATCATCTGTGTGGACACTGAAGTTTCTGTCACATACTACTGTGCTCACAGACTACAG",
            "read frame": 0,
            "strand": True
        },
        "start": 666355,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "3609.9.164*01",
        "confirmed": True,
        "end": 602739,
        "family": "3609v",
        "framed": True,
        "functionality": "O",
        "gene": "3609.9.164",
        "id": "3609.9.164*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-9*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-9",
        "length": 301,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTACTCTGAAAGAGTCTGGCCCTGGTATATTGCAGCCCTCCCAGACCCTCAGTCTGACCTGTTCTTTCTCTGTGTTTTCACTGAGCACTTTTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTGGAGTGGCTGGCACACATTTATTGGGATGAGGACAAGCACTATAAACCATCCTTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCACTGTGGACACTGCAGATACTGCCACATACTACTCTGCTCGAAGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 602438,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.15.75*01",
        "confirmed": True,
        "end": 1738313,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "PG.15.75",
        "id": "PG.15.75*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-1*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-1",
        "length": 291,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTTAGCTTCTAGAATCTGGAGGCTCTTTGATATAGCCTGGAGGATTCTTTAAACTCTCTTCTAAAGCTGCTGGATTCACCTTCAATGATTACTGGATGCACTGGTTTTGATAAGGCTCAGGGAATCGGCTAGAGTGGTTAGCAGATATAAAATATGTGAAAAGTGTAAAAAACCATGCAGAGTCTGTGAAGGGAAGATTTGCCACTACAGAGACAATTCTAAAAACTTATTTTATCTACAAATGAACAGTCTGAGAAATGAGGACATTGCCCCTTATTACTCTATGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1738022,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.10pg.16*01",
        "confirmed": True,
        "end": 2401583,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.10pg.16",
        "id": "7183.10pg.16*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-10*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-10",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTAAAGCTGGTGGAGTCTAGGGGAGGCTTAGTGTAGCCTGGAAGGTCCGTGATACGCTCATGTGCAGCCTCTGGATTCACTGACAGTGACTAATAGTTGGCCTAGGTTTGCCAAGTTCCAAAGAAGGAGCTGGAATGGGGGGCATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGAAGAAGTGGGTTGGACATAACATATTCAGAATCAATATTTAGAGATTCTAATCCTTGAAGACATCACTTTTGACCAAGTATATATGCAACATGTTACTTAGGTT",
            "read frame": 0,
            "strand": True
        },
        "start": 2401287,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.11pg.19*01",
        "confirmed": True,
        "end": 2384027,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.11pg.19",
        "id": "7183.11pg.19*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-11*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-11",
        "length": 265,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTAGGGTCAGGGCAGCCTGGAGGGTCCCTGAAACACTCCTGTGCAGCCTCTGTAGTCACTGTGAGTGACTACTGAATGACCTGGGTCCTTCAGGCTCTAAAGAAGGGGCTGGAGAGGGTGGAAATAATTTTTAATGGTGGAGGTAGCACCTATTATCCAGACACCATGAAGGGCTGATTCACCATCTACAGAGATGATGCCAGAAACACACTTTACCTGAAAATAAACAGTCTGAGGTCTGAGTACACAGCCA",
            "read frame": 0,
            "strand": True
        },
        "start": 2383762,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.12.20*01",
        "confirmed": True,
        "end": 2369027,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.12.20",
        "id": "7183.12.20*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-12*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-12",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGAAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTGACTATTACATGTATTGGGTTCGCCAGACTCCAGAGAAGAGGCTGGAGTGGGTCGCATACATTAGTAATGGTGGTGGTAGCACCTATTATCCAGACACTGTAAAGGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCCGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2368731,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.13pg.24*01",
        "end": 2346396,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.13pg.24",
        "id": "7183.13pg.24*01-C57BL/6",
        "imgt family name": "IGHV5",
        "length": 295,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGTTCAGGGTGAGGTAGAGCTGCTGAAATCCAAGGAAGGCATAGTGATACCTAAAAGCCCCTGAAACTCTCATGATAAACCTCTGGATTCACATTCAGTGATTACTATATAAGCAGGGTCCTCCACCCTCTAGGGAAGGGCCTTCAGTGGTGAGCATACATTAATAGAGACAGCTCCATCAACTGTGCTGAAGCCATGAAAAGCCAATTCACCATCACCAGAGACAATGCCAAGAACAACCTCTACCTGGATATTTGCAGACTGAGGACACAGCCATGTTTTCTTGTGCAGGATA",
            "read frame": 0,
            "strand": True
        },
        "start": 2346101,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.14.25*01",
        "confirmed": True,
        "end": 2335039,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.14.25",
        "id": "7183.14.25*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-9-1*02",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-9-1",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GACGTGAAGCTGGTGGAGTCTGGGGAAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGCCATGTCTTGGGTTCGCCAGACTCCAGAGAAGAGGCTGGAGTGGGTCGCATACATTAGTAGTGGTGGTGATTACATCTACTATGCAGACACTGTGAAGGGCCGATTCACCATCTCCAGAGACAATGCCAGGAACACCCTGTACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTACAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2334745,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.15pg.26*01",
        "end": 2325722,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.15pg.26",
        "id": "7183.15pg.26*01-C57BL/6",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-7-6",
        "length": 298,
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGCCCATTGTGAGGTGAAGCTGGTGGAGTCTAAGGGAGGTTTAGTGTAGCCTGGAAGGTCCATGATACTCTACTGTGCAGCCTCTGGATTCACTGTCAGTGACTACTGGTTGGCTGGGTTTGCCAGGCTCCAAAGAAGGAGCTGGAATGGGTGGGGACATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGTAAAGTGGGTTGGACATGACATATTCAGAATCAATATTTAGAGATTCTAATCCTTGAAGACATCACTTTCGACCAAGTATATATGCAAGATGTT",
            "read frame": 0,
            "strand": True
        },
        "start": 2325424,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.16.27*01",
        "confirmed": True,
        "end": 2308897,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.16.27",
        "id": "7183.16.27*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-12-4*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-12-4",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGAAGCTCGTGGAGTCTGGGGGAGGTTTAGTGAAGCCTGGACAGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAACTATTACATGTCTTGGGTTCACCAGACTCCGGAGAAGAGGCTGGAGTGGGTTGCATACATTAGTAGTAGTGGTGTTAGCACCTATTATCCAGACAATGTAAAGGGCCGATTCGCCATCTCCAGAGACAATGCCAAGAACACCCTTTACCTGCAAATGACCAGTCTGAAGTCAGAGGACACGGCCTTGTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2308603,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.17pg.31*01",
        "end": 2289766,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.17pg.31",
        "id": "7183.17pg.31*01-C57BL/6",
        "imgt family name": "IGHV5",
        "length": 294,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGTTCAGGGTGAGGTAGAGCTGCTGAAATCCAAGGAAGGCATAGTGATGCCTAAAAGTCCCTGAAACACTCATGTTAAACCTCTGGATTCACATTCAGTGATACTATATAAGCAGGGTCCTCCACCCTCTAGGGAAGGGTCTTCAGTTGTGAGCATACATTAATAGAGACAGCTCCATCAACTGTGCTGAAGCCATGAAAAGCCAACTCACCTTCACCAGAGACAATGCCAAGAACAACCTCTACCTGGAAATTTGCAGGCTGAGGACACAGCCATGTTTTCTTGTGCAGGATA",
            "read frame": 0,
            "strand": True
        },
        "start": 2289472,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.18.35*01",
        "confirmed": True,
        "end": 2244503,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.18.35",
        "id": "7183.18.35*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-15*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-15",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTGACTACGGAATGGCGTGGGTTCGACAGGCTCCAAGGAAGGGGCCTGAGTGGGTAGCATTCATTAGTAATTTGGCATATAGTATCTACTATGCAGACACTGTGACGGGCCGATTCACCATCTCTAGAGAGAATGCCAAGAACACCCTGTACCTGGAAATGAGCAGTCTGAGGTCTGAGGACACGGCCATGTATTACTGTGCAAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2244207,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.19.36*01",
        "confirmed": True,
        "end": 2232623,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.19.36",
        "id": "7183.19.36*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-16*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-16",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGAAGCTGGTGGAGTCTGAGGGAGGCTTAGTGCAGCCTGGAAGTTCCATGAAACTCTCCTGCACAGCCTCTGGATTCACTTTCAGTGACTATTACATGGCTTGGGTCCGCCAGGTTCCAGAAAAGGGTCTAGAATGGGTTGCAAACATTAATTATGATGGTAGTAGCACCTACTATCTGGACTCCTTGAAGAGCCGTTTCATCATCTCGAGAGACAATGCAAAGAACATTCTATACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCACGTATTACTGTGCAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2232327,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "AF120460",
        "allele": "7183.1b*01",
        "confirmed": True,
        "end": 297,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.1b",
        "id": "7183.1b*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-2*03",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-2",
        "length": 297,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGAAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGAGTCCCTGAAACTCTCCTGTGAATCCAATGAATACGAATTCCCTTCCCATGACATGTCTTGGGTCCGCAAGACTCCGGAGAAGAGGCTGGAGTTGGTCGCAGCCATTAATAGTGATGGTGGTAGCACCTACTATCCAGACACCATGGAGAGACGATTCATCATCTCCAGAGACAATACCAAGAAGACCCTGTACCTGCAAATGAGCAGTCTGAGGTCTGAGGACACAGCCATGTATTACTGTGCAAGACGA",
            "read frame": 0,
            "strand": True
        },
        "start": 0,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.1pg.1*01",
        "confirmed": True,
        "end": 2498329,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.1pg.1",
        "id": "7183.1pg.1*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-1*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-1",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGGGTCCCGGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGCCATGTCTTGGGTCCGCCAGACTCCGGAGAAGAGGCTGGAGTGGGTCGCAGCCATTAGTACTGATGGTAGTTTCATCTACTAACCAGACACTGTAAAAGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTTTCTGCAAATGAGCAGTCTAAGGTATGAGGACACGGCCATGTATTACTGTTTGAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2498033,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.2.3*01",
        "confirmed": True,
        "end": 2492752,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.2.3",
        "id": "7183.2.3*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-2*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-2",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGAGTCCCTGAAACTCTCCTGTGAATCCAATGAATACGAATTCCCTTCCCATGACATGTCTTGGGTCCGCAAGACTCCGGAGAAGAGGCTGGAGTTGGTCGCAGCCATTAATAGTGATGGTGGTAGCACCTACTATCCAGACACCATGGAGAGACGATTCATCATCTCCAGAGACAATACCAAGAAGACCCTGTACCTGCAAATGAGCAGTCTGAGGTCTGAGGACACAGCCTTGTATTACTGTGCAAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2492456,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.20.37*01",
        "confirmed": True,
        "end": 2212000,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.20.37",
        "id": "7183.20.37*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-17*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-17",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTGACTATGGAATGCACTGGGTTCGTCAGGCTCCAGAGAAGGGGCTGGAGTGGGTTGCATACATTAGTAGTGGCAGTAGTACCATCTACTATGCAGACACAGTGAAGGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTTCCTGCAAATGACCAGTCTGAGGTCTGAGGACACGGCCATGTATTACTGTGCAAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 2211706,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.21pg.38*01",
        "confirmed": True,
        "end": 2195382,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.21pg.38",
        "id": "7183.21pg.38*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-18*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-18",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GATGTAGCTGGTGGAGTCTAGGAGAGGCTTAGTACAGCCTGGAAGGTCCCTGAACTCTCATGAGCAGCCTCTGGATTCACTTTTAGTAATTATGTCATGGCCTGGGTCCAACAAGCTCCAGAGAAGGGGCTGGAGTGGGTTGCATACATTAGTAGTGGCAGTAGTACCTTCTATTATGCAGACACAGTTAAGGGCCCATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTTCCAGCAAATGAGCAGTCTAAGGTCTGAGGACAGAGACAAAAAGATAAAAATATATGGCCA",
            "read frame": 0,
            "strand": True
        },
        "start": 2195086,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.3pg.5*01",
        "confirmed": True,
        "end": 2481289,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.3pg.5",
        "id": "7183.3pg.5*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-3*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-3",
        "length": 265,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGATAGGGTCAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGAGTCACTGTGAGTGACTACTGAATGACCTGGGTCCTTCAGGCTCCAAAGAAGGGGCTGGAGAGGGTGGCAATAATTTTTAATGGTGGAGGTAGCACCTATTATCCAGAAACCATGAAGGGCTGATTCACTATCTACAGAGATGATGCCACAAACACACTTTACCTGAAAATAAACAGTCTGAGGTCTTAGTACACAGCCC",
            "read frame": 0,
            "strand": True
        },
        "start": 2481024,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.4.6*01",
        "confirmed": True,
        "end": 2473810,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.4.6",
        "id": "7183.4.6*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-4*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-4",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGCCATGTCTTGGGTTCGCCAGACTCCGGAAAAGAGGCTGGAGTGGGTCGCAACCATTAGTGATGGTGGTAGTTACACCTACTATCCAGACAATGTAAAGGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACAACCTGTACCTGCAAATGAGCCATCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2473514,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.6pg.9*01",
        "confirmed": True,
        "end": 2458263,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.6pg.9",
        "id": "7183.6pg.9*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-5*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-5",
        "length": 268,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTGGAGTCTGTGGGAGGGTTAGAGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGAGTCACTGTCAGTGAGTACTGAATGACCTGGGTCCTTCAGGCTCCAAAGAAGGAGCTGGAAAGGGTGGCAATAATTTTTAATGGTGGAGGTAGCACCTATTATCCAGAAACCATGAAGGGCCAATTCACCATCTGCAGAGATGATGCCAAAAACATACTTTACCTAGAAATAAACACTCTGTGGTCTGAGT",
            "read frame": 0,
            "strand": True
        },
        "start": 2457995,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.7.10*01",
        "confirmed": True,
        "end": 2445750,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.7.10",
        "id": "7183.7.10*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-6*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-6",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGCAGCTGGTGGAGTCTGGGGGAGACTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGGCATGTCTTGGGTTCGCCAGACTCCAGACAAGAGGCTGGAGTGGGTCGCAACCATTAGTAGTGGTGGTAGTTACACCTACTATCCAGACAGTGTGAAGGGGCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2445454,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.8pg.14*01",
        "confirmed": True,
        "end": 2416289,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.8pg.14",
        "id": "7183.8pg.14*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-8*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-8",
        "length": 277,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTGGAGTCTGTGGGAGGCTTAGTGCAGCCTGGAGGGTCGCTGAAACTCTCCTGGGCAGCCTCTGGATTAACTCACTGACAACTGAATGACCTGGGTCCTTCAGGCTCCAAGGAAGGGGCTGGAGAGGGTGGCAATAATTTTTAGTGGTGGAGGTAGCACCTACTATATAGACGCAATGAAGGGCCGATTCACCGTCTCCAGAGATGATAACAAAAACACACTTTACCTGAAAATAAACAGTCTGAGGTCTGAGTACACAGCCATG",
            "read frame": 0,
            "strand": True
        },
        "start": 2416012,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.9.15*01",
        "confirmed": True,
        "end": 2409487,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.9.15",
        "id": "7183.9.15*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-9*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-9",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGATGCTGGTGGAGTCTGGGGGAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATACCATGTCTTGGGTTCGCCAGACTCCGGAGAAGAGGCTGGAGTGGGTCGCAACCATTAGTGGTGGTGGTGGTAACACCTACTATCCAGACAGTGTGAAGGGTCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCAGTCTGAGGTCTGAGGACACGGCCTTGTATTACTGTGCAAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2409191,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.1.11*01",
        "confirmed": True,
        "end": 2436570,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "PG.1.11",
        "id": "PG.1.11*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-7*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-7",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTGGAGTCTAGGGGAGGCTTAGTGTAGCCTGGAAGGTCCGTGATATGCTCATGTGCAGCCTCTGGATTCACTGACAGTGACTAATAGTTGGCCTAGGTTTGCCAAGTTCCAAAGAAGGAGCTGGAATGGGGGGCATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGAAGAAGTGGGTTGGACATAACATATTCAGAATCAATATTTAGAGATTCTAATCCTTGAAGACATCACTTTTGACCAAGTATATATGCAACATGTTACTTAGGTT",
            "read frame": 0,
            "strand": True
        },
        "start": 2436274,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.14.73*01",
        "confirmed": True,
        "end": 1751080,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "PG.14.73",
        "id": "PG.14.73*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-21*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-21",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAGGCTGGTGGAATCTGGAGGCAGCTTGATAGAGCATGAAGGGTACATTCAACTCTTTTGTCAAGCTTCTGGATTCACCCTCAGTGGTTACTGGATGCACTGGATTTGCCAAGCCCCAGGGAAGGGGCCAGAGTGGTTAGCAAATATAAAATATGATGGGAGTGAAAAATACTATGCAGTGTCTATGAAGGGGTGATTTGCCATCTCCAGAGACCTTCCTAAGAACTTTCTTTATCTGCAAATGAGCAATTTGAGAAATGAGGACACTGCAATGTATTACTGTGCAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1750784,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.2.21*01",
        "end": 2356544,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "PG.2.21",
        "id": "PG.2.21*01-C57BL/6",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-13",
        "length": 289,
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGAGGTGAAGCTGGTGAAGTCTAAGGGGAGGCATAGTGCAGCCTAGAAGGTCCATGATACTCTACTGTGCAGCCTCGGATTCACTGTAAGTGACGACTGGTTTGTCCGTGTTTGCCAGGCTCCAAAGAAGGGGCTGCAGTGGGGGATGGGAATAATTTTTCATGGTTGTGGTAGCCCCTCTTATGCAGACACCTTGAAGAAGTGGGTTGGACGTAACATATTCAGAATCAATATTTAAAGATTCTAATCCTTGAAGATATCACTTTTGACCAAGTATATATGAACCAT",
            "read frame": 0,
            "strand": True
        },
        "start": 2356255,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.7.41*01",
        "confirmed": True,
        "end": 2186386,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "PG.7.41",
        "id": "PG.7.41*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-19*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-19",
        "length": 297,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATATCCAGTGTGAGGTGCAGCTGGTGGGGTCCAGGAGAGTCTTAATGCAGCCTAGAAAAACCTTGAAACTCCCCTCTGCAGCACATGGAATCACCTTCAGCGACTGACATGGTACTCCGGGCTCCAGGGAAGGGGCTGGAGTGGGATACATTCATTAATAGTAAGGGTAGTTACATCAACTATGCCAAAGCTATAAAAGAACAATTCACTATCTCCAGAGACAATACCAAGAACACTTTGTACCTGGAAATTAGCAGTTTGAAGTCTGAAAACACAGTCATTTATTACTGTACAACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2186089,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.1.85*01",
        "confirmed": True,
        "end": 1602034,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.1.85",
        "id": "J558.1.85*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-2*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-2",
        "length": 297,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CCTGAGCTGGCAAGGCCTTGGGCTTCAGTGAAGATATCCTGCCAGGCTTTCTACACCTTTTCCAGAAGGGTGCACTTTGCCATTAGGGATACCAACTACTGGATGCAGTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATCGGGGCTATTTATCCTGGAAATGGTGATACTAGTTACAATCAGAAGTTCAAGGGCAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAACTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCATGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1601737,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.10pg.100*01",
        "confirmed": True,
        "end": 1440428,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.10pg.100",
        "id": "J558.10pg.100*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-13*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-13",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTTCAGCAGTCTGGACCAGAGCTGGTAATACCTGGGGCTTAGGTGAAGTTGTCCTGCAAGGCTTCTGGCTACAATTTTAATGACTATGAAATTCAATGGGTGAAGCAGAGTCTGAAGCAGGGACTGGAATGGATTGGAGCTATTCATCCTGAAAATGGTGGTATTACCTACAATCAGAAGTTCAAAGGCAAGGCCACATTTACTGTAGACACATCCTCCAACACAGCCTACATGCAACTCAGAAGCCTGACATCTGAGGACACTGCTGACTATTATTGTGAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1440134,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.11pg.101*01",
        "confirmed": True,
        "end": 1424536,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.11pg.101",
        "id": "J558.11pg.101*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-14*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-14",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGTTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACATTCACTAGCTATGTTATGCACTGGGTGAAGCAGAAGCCTGGGCAGGGCCTTGAGTGGATTGGATATATTTATCCTTACAATGATGGTACTAAGTACAATGAGAAGTTCAAAGGCAAGGCCACACTGACTTCAGACAAATCCTCCAGCACAGCCTACATGGAGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1424242,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.12.102*01",
        "confirmed": True,
        "end": 1413755,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.12.102",
        "id": "J558.12.102*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-15*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-15",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAACTGCAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGACGCTGTCCTGCAAGGCTTCGGGCTACACATTTACTGACTATGAAATGCACTGGGTGAAGCAGACACCTGTGCATGGCCTGGAATGGATTGGAGCTATTGATCCTGAAACTGGTGGTACTGCCTACAATCAGAAGTTCAAGGGCAAGGCCATACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCCGTCTATTACTGTACAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1413461,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.13.103*01",
        "confirmed": True,
        "end": 1405236,
        "family": "J558",
        "framed": True,
        "functionality": "O",
        "gene": "J558.13.103",
        "id": "J558.13.103*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-16*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-16",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAGCAGTCTGGACCTAAGGTAGTGAATGCTGGGGCTTCCGTGAAGCTGTCCTGCAAGTCTTCTGGTTACTCATTCAGTAGATACAAAATGGAATGTGTGAAACAGAGCCATGTAAAGAGCCTTGAGTGGATTGAACATATTAATCTTTTCAATGGTATTACTAACTACAATGGAAACTTTAAAAGCAAGGCCACATTGACTGTAGACATATCCTCTAGCACAGCCTATATGGAGCTTAGCAGATTGACATCTGAAGACTCAGAGGTATATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1404942,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.14pg.104*01",
        "confirmed": True,
        "end": 1399140,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.14pg.104",
        "id": "J558.14pg.104*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-17-1*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-17-1",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTCCACCTGCAGCAGTCTCTACCTAAGGTAGTGAAGGCTGGGCCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGATCCTTCAAAGGATTGAATATGTTAATCCTTATAATGGTGGTACTGGCTACAATGAAAAGTTCAAGGACAAGGCCACATTGACTGCAGACAAATCCTTCAGCACAGCCTATATGCAATTCAGCAGCCTGACATCTGAGGACTCTCTGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1398846,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.15pg.105*01",
        "confirmed": True,
        "end": 1394359,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.15pg.105",
        "id": "J558.15pg.105*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-17*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-17",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTATTCATTCAGTGACTACTACATGGAATGGGTGAAGCAGAGCCATAGAAAGAGCCTTGAATGTATTGGAGAAATTAATCCTTACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCTTCCAGCACAGCGTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTTCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1394066,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.16.106*01",
        "confirmed": True,
        "end": 1388476,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.16.106",
        "id": "J558.16.106*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-18*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-18",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATACCCTGCAAGGCTTCTGGATACACATTCACTGACTACAACATGGACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGATATTAATCCTAACAATGGTGGTACTATCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACACTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1388182,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.17pg.107*01",
        "confirmed": True,
        "end": 1366933,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.17pg.107",
        "id": "J558.17pg.107*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-19-1*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-19-1",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGAAGCAGCCTGGAACTGTGGTGGTGAAACCTGGGGCTTCAGTGAAGATATCCTGCCAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTGTAGTGGATTGGACTTATTATTCCTTACAATGGTGATACTGGCTACAACCAGAAGTTCAAGGGCAAAGCAACATTGACTGTAGACAAGTCCTCCAGCACAGCCAACATGGAGCTCTGCAGCCTGACATCTGAGGACTCTACAGTCTATTACCGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1366639,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.18.108*01",
        "confirmed": True,
        "end": 1362432,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.18.108",
        "id": "J558.18.108*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-19*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-19",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAACAGTCTGGACCTGTGCTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGTAAGGCTTCTGGATACACATTCACTGACTACTATATGAACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGTTATTAATCCTTACAACGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTTGACAAGTCCTCCAGCACAGCCTACATGGAGCTCAACAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1362138,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.19.109*01",
        "confirmed": True,
        "end": 1347308,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.19.109",
        "id": "J558.19.109*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-20*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-20",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGATTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTTACTGGCTACTTTATGAACTGGGTGATGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACGTATTAATCCTTACAATGGTGATACTTTCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCTAGCACAGCCCACATGGAGCTCCGGAGCCTGACATCTGAGGACTCTGCAGTCTATTATTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1347014,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.2.88*01",
        "confirmed": True,
        "end": 1583973,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.2.88",
        "id": "J558.2.88*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-4*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-4",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGGGCTGAACTGGCAAGACCTGGTGCCTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTTACTAGCTACACGATGCACTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATTGGATACATTAATCCTAGCAGTGGTTATACTAAGTACAATCAGAAGTTCAAGGACAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAACTGAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1583679,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.21pg.111*01",
        "confirmed": True,
        "end": 1330672,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.21pg.111",
        "id": "J558.21pg.111*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-21*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-21",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCATGGCTTCTGGTTATTCATTCAGTGACTACTACATGCACTGAGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGATATATTAACCCTAACAATGGTTGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCGTACATGGAGCTCCACAGCCTGACATCTGAGGACTCTTCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1330378,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.22.112*01",
        "confirmed": True,
        "end": 1324807,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.22.112",
        "id": "J558.22.112*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-22*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-22",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACATTCACTGACTACAACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGATATATTAACCCTAACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAAACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCGGAGGATTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1324513,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.23.113*01",
        "confirmed": True,
        "end": 1306630,
        "family": "J558",
        "framed": True,
        "functionality": "O",
        "gene": "J558.23.113",
        "id": "J558.23.113*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-23*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-23",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAACTGCAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCGGGCTACACATTTACTGACTATGAAATGCACTGTGTGAAGCAGACACCTGTGCACGGCCTGGAATGGATTGGAGCTATTGATCCTGAAACTTGTGGTACTGCCTACAATCAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCCGTCTATTACTGTACAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1306336,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.24pg.114*01",
        "confirmed": True,
        "end": 1298152,
        "family": "J558",
        "framed": True,
        "functionality": "O",
        "gene": "J558.24pg.114",
        "id": "J558.24pg.114*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-24*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-24",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCGGCTGCAGCAGTCTGGACCTAAGGTAGTGAATGCTGGGGCTTCCGTGAAGCTGTCCTGCAAGTCTTCTGGTTACTCATTCAGTAGATACAAAATGGAATGTGTGAAACAGAGCCATGGAAAGAGCCTTGAGTGGATTGAACATATTAATCTTTTCAAAGGTGTTACTAACTACAATGGAAAGTTCAAAAGCAAGGCCACATTGACTGTAGACATATCCTCTAGCACAGCCTATATGGAGCTTAGCAGATTGACATCTGAAGACTCAGAGGTATATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1297858,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.25pg.115*01",
        "confirmed": True,
        "end": 1289671,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.25pg.115",
        "id": "J558.25pg.115*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-25*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-25",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCGTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACATTATGAACTTGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTTGAGAAATTAATCCTTACAATGGTGGTACTAACTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGAGTACATGGAGCTCCACAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1289377,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.26.116*01",
        "confirmed": True,
        "end": 1282708,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.26.116",
        "id": "J558.26.116*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-26*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-26",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAACAATCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGTAAGGCTTCTGGATACACGTTCACTGACTACTACATGAACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGATATTAATCCTAACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1282414,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.27pg.117*01",
        "confirmed": True,
        "end": 1266464,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.27pg.117",
        "id": "J558.27pg.117*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-27*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-27",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "AAGGTCAAGCTGTAGCAGTCTGGACCTGAGCTGGTGAAGTCTGGGGCTTCAGAGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGATCATTGAGTAGATTGGACTTATTATTCCCTACAATGGTGATACTGGCTACAACCAGAAGTTTAAGGGCAAGGCCACATTGACTGTAGACTAGTCCTTCAGCACAGCCTACATGGAGCTCCACAGCCTGAAATCTTAGGACTCTGTGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1266170,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.28pg.118*01",
        "confirmed": True,
        "end": 1261659,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.28pg.118",
        "id": "J558.28pg.118*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-28*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-28",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGATCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGTGCTTCAGTGAAGATATCCTTCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGAACTGGATGAAGTAGAGCCACGGAAATAGCCTTGAGTGGATTGGATATATTGATCCTTACAATGGTGGTACTAGCTACAACCAGAAATTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCTTCCAGCACAGCCTACATGCATCTCAACAACCTGACATCTGAGGACTCTGCAGTCTATTACTTTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1261365,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.29pg.119*01",
        "confirmed": True,
        "end": 1258759,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.29pg.119",
        "id": "J558.29pg.119*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-29*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-29",
        "length": 285,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCACCTGTAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTCCTGGTTACTCATTCATTGGCTACTACATGCACTGAGTAAAGTCCTGTAAAAAGCCTTTAAAGGATTGGATATATTAATCCAGTGGTGGTACTGGCTACAATGAAAAGTTCAAGGGCAAGGCCACATTGACTGCAGACAAATCCTTCAGCACAGCCTATATGCAATTCAGCAGCCTGACATCTGAGGACTCTGTGGACTCTTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1258474,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.3.90*01",
        "confirmed": True,
        "end": 1557779,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.3.90",
        "id": "J558.3.90*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-5*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-5",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTCAGCTCCAGCAGTCTGGGACTGTGCTGGCAAGGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGACTTCTGGCTACACATTTACCAGCTACTGGATGCACTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATAGGGGCTATTTATCCTGGAAATAGTGATACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCAAACTGACTGCAGTCACATCCGCCAGCACTGCCTACATGGAGCTCAGCAGCCTGACAAATGAGGACTCTGCGGTCTATTACTGTACAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1557485,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.30pg.120*01",
        "confirmed": True,
        "end": 1253981,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.30pg.120",
        "id": "J558.30pg.120*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-30*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-30",
        "length": 292,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGTTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGACTACTACATGCACTGGGTGAAGCAAAGTCCTGAAAAGAGTCTTGAGTGGATTGGAGAGATCAATCCTAGCACTGTTGTTACTAACTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAAGAGCCTGACATCTGAAGACTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1253689,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.31.121*01",
        "confirmed": True,
        "end": 1241816,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.31.121",
        "id": "J558.31.121*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-31*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-31",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAATATCCTCGATTGGATTGGATATATTTATCCTTACAATGGTGTTTCTAGCTACAACCAGAAATTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCTAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1241522,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.32pg.122*01",
        "confirmed": True,
        "end": 1235857,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.32pg.122",
        "id": "J558.32pg.122*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-32*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-32",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCAGCAAGACTTCTGGTTAAACTTTTACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTTGATATATTGATCCTTACGATGGTGTTACTAGCTACAATAAGAAGTTCAAGAGAAAGGCCACGTTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCTGAAACCTGACATCTGAGGAAACTGAAGTCTATTACTGTGTAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1235563,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.33pg.123*01",
        "confirmed": True,
        "end": 1227429,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.33pg.123",
        "id": "J558.33pg.123*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-33*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-33",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAACAATCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGTAAGGCTTCTGGATACACATTCACTGATTACTACATGCACTGGGTGAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGATATTAATCCTTACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAACAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1227136,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.34.124*01",
        "confirmed": True,
        "end": 1219890,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.34.124",
        "id": "J558.34.124*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-34*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-34",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAACAGTCTGGACCTGAGTTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACATTCACTGACTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGATATATTTATCCTAACAATGGTGGTAATGGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1219596,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.35pg.125*01",
        "confirmed": True,
        "end": 1195776,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.35pg.125",
        "id": "J558.35pg.125*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-35*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-35",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGAAGCAGCCTGGAACTGTGGTGGTGAAACCTGGGGCTTCAGTGAAGATATCCTGCCAGGCTTCTGGTTACTCATTTACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGAAAAGAGCCTTTAGTGGATTCTACTTATTATTCCTTACAATGGTGATACTAGCAACAACCAGAAGTTCAAGGGCAAAGCAACATTGACTGTAGACAAGTCCTCCAGCACAGCCAACATGGAGCTCTGCAGCCTGACATCTGAGGACTCTACGGTCTATTACCGTGCAAGACT",
            "read frame": 0,
            "strand": True
        },
        "start": 1195480,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.36.126*01",
        "confirmed": True,
        "end": 1191189,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.36.126",
        "id": "J558.36.126*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-36*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-36",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAACAGTCTGGACCTGTGCTGGTGAAGCCTGGGCCTTCAGTGAAGATATCCTGTAAGGCTTCTGGATTCACATTCACTGACTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACTTGTTTATCCTTACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGGAGCTAAACAGCCTGACTTCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1190895,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.37.127*01",
        "confirmed": True,
        "end": 1174839,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.37.127",
        "id": "J558.37.127*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-37*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-37",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTTACTGGCTACTTTATGAACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACGTATTAATCCTTACAATGGTGATACTTTCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCTAGCACAGCCCACATGGAGCTCCTGAGCCTGACATCTGAGGACTTTGCAGTCTATTATTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1174545,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.38pg.128*01",
        "confirmed": True,
        "end": 1161066,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.38pg.128",
        "id": "J558.38pg.128*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-38*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-38",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTCGTGAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTAGCTACTAAATGCACTGGGTGAAGCAAAGCCATGGAAAGAGCCTTGAGTGGATTGGACTTATTATTCCTTACAATGGTGATACTGGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACTAGTCCTTCAGCACAACCTACATGGAGCTCCGCAGCCTGAAATCTGAGGACTCTGCGGTCTATTACTCTACAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1160773,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.39.129*01",
        "confirmed": True,
        "end": 1156478,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.39.129",
        "id": "J558.39.129*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-39*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-39",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGTTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGCGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGACTACAACATGAACTGGGTGAAGCAGAGCAATGGAAAGAGCCTTGAGTGGATTGGAGTAATTAATCCTAACTATGGTACTACTAGCTACAATCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACCAATCTTCCAGCACAGCCTACATGCAGCTCAACAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1156184,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.4.93*01",
        "confirmed": True,
        "end": 1532613,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.4.93",
        "id": "J558.4.93*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-7*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-7",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGGGCTGAACTGGCAAAACCTGGGGCCTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTTACTAGCTACTGGATGCACTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATTGGATACATTAATCCTAGCAGTGGTTATACTAAGTACAATCAGAAGTTCAAGGACAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTGAGCAGCCTGACATATGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1532319,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.40pg.130*01",
        "confirmed": True,
        "end": 1146136,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.40pg.130",
        "id": "J558.40pg.130*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-40*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-40",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAACTGCAACAGTCTGGATCTAAGGTAGTGAATGCTGGGGCTTCCATGAAGCTGTCCTGCAAATATTCTGGTTACTCATTTAGTAGATACAAAATGGAATGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTCGATTGGAGATATTAATCTTTCCAATGGTGGTACTAACTACAATGGAAAGTTCAAAAGCAAGGCCACATTGACTGTAGATATATCCTCTAGCACAGCCTATATGGAGCTAGCATATTGACATCTGAGGTCTCTGCAGTCTCTCACCATGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1145843,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.41pg.131*01",
        "confirmed": True,
        "end": 1138706,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.41pg.131",
        "id": "J558.41pg.131*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-41*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-41",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TCTCCACATAGGTGTCCACTCTGAGGTCCACCTGAAGCAGTCTGGACCTAAGGTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGTCCTGTAAAAAGCCTTCAAAGGATTGGATATATTATTCCTTACAGTGGTGGTACTGGGTACAATGAAAAGTTCAAGGGCAAGGCCACATTGTCTGCAGACAAATTCTTCAGCACAACCTATATGTACTTCAGCAGCCTGACATCTGAGGACTCTGTGGACTCTTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1138412,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.42.132*01",
        "confirmed": True,
        "end": 1133964,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.42.132",
        "id": "J558.42.132*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-42*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-42",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGAACTGGGTGAAGCAAAGTCCTGAAAAGAGCCTTGAGTGGATTGGAGAGATTAATCCTAGCACTGGTGGTACTACCTACAACCAGAAGTTCAAGGCCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAAGAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1133670,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.43.133*01",
        "confirmed": True,
        "end": 1125127,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.43.133",
        "id": "J558.43.133*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-43*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-43",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCAAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAAAGTTCTGAAAAGAGCCTTGAGTGGATTGGAGAGATTAATCCTAGCACTGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTAACTGTAGACAAGTCATCCAGCACAGCCTACATGCAGCTCAAGAGCCTGACATCTGAGGACTCTGCTGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1124833,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.44pg.134*01",
        "confirmed": True,
        "end": 1107697,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.44pg.134",
        "id": "J558.44pg.134*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-44*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-44",
        "length": 295,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGGGCTGAAATGGTGAGGCCTGGGGCCTCAGTGAAGATGTCCTGCAAGGCTTTTGGCTACACATTCACTGACTATGGCATGCATTGGGTGAAGCATAGTCACGGAAGGCATCCTTGAGAGGATTGGAAATATTAATACTTACTATGACAATACTGGCTACAATGAAAAGTTCAAGGGCAAGTCCAAATTGACTGTAGACAAATCCTCCAGCACAGCCTATGTGGAGTTTAGCAGAATGACATCTGAGGATTCTGTAGTCTATTACTGTGAAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1107402,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.45pg.135*01",
        "confirmed": True,
        "end": 1101022,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.45pg.135",
        "id": "J558.45pg.135*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-45*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-45",
        "length": 295,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAGATGCAGCAGTCTGGGGGATTAGGGGATGAAGCCTGTGTTCTCAGTGAAGATGTCCTGTAAGGGTTCTGGCTACAGCTTTACCAACTACTATATGCCCTGGGTAAAATAGTGGACTGGACACAACCTTGAGTGGATTGGATGGATTCATCCTGGAAATGGTGATACTTACTACAATCAAAAGTTCAAGGGAAAGGCAACACTGACCAAGTACAAATTCTCCAGCACAGCCTACTTACATCACAACAGCCTGACATCTGAGCACCCAGTAGTTTATAAATATGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1100727,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.46pg.136*01",
        "confirmed": True,
        "end": 1093384,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.46pg.136",
        "id": "J558.46pg.136*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-46*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-46",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGGTGCAGCTGTCTGCAGCTGAGCTGGTGAAGCCTGGGAGTCCAGTGAAGCTGTCCTGCAAAGCTTCTGGCTACACCGTCAATGACAACTATATGGAGCAGGTAAAGCAGAGGCCTGGACAGAGCATGGAATGGATTGGATAGATTCATTTTGTATATGGTGGTACTTAATACAATGAAAAGTTCTAGGGCAAGTCCACATTAACTGTAGAAAAATCCTCCAACACAGCCTACATGGAACTCAACAGCTCGACATCTGAGGACTCTGTAGTTTATTACTGTGCATGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1093090,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.47.137*01",
        "confirmed": True,
        "end": 1079969,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.47.137",
        "id": "J558.47.137*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-47*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-47",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAGCTGCAGCAGTCTGGGGCTGAGCTAGTGAAGCCTGGAGCCTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTCACTACCTATCCTATAGAGTGGATGAAGCAGAATCATGGAAAGAGCCTAGAGTGGATTGGAAATTTTCATCCTTACAATGATGATACTAAGTACAATGAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGAAAAATCCTCTAGCACAGTCTACTTGGAGCTCAGCCGATTAACATCTGATGACTCTGCTGTTTATTACTGTGCAAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1079675,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.48pg.140*01",
        "confirmed": True,
        "end": 1032850,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.48pg.140",
        "id": "J558.48pg.140*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-48*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-48",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TAGGTCCAATGGCAGGAGTCAGGGACTGAGCTGGTGAGATCTGGGGCCTCAGTAATGATGTCCTGCAAGGCTTCTGGATACACATTCAGTAACTACACTATGCACTGGGTAAAGCAGAGTCATGGAAAAGGCCATGAAAGGATTGGATATATTGATAATTACTATTGTAGCACTGACTACAGTGAAAAGTTCAAGATCAAGGCCACATTGACTGTAAACAAATCCTGCAGAACAGCCTATGTCAAGCTCAGCAGACTGACATCTGAGGACTCTGCAGTCTATTATTGTGTAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1032556,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.49.141*01",
        "confirmed": True,
        "end": 1015854,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.49.141",
        "id": "J558.49.141*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-49*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-49",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGCGTGAGCTGCAGCAGTCTGGAGCTGAGTTGGTGAGACCTGGGTCCTCAGTGAAGTTGTCCTGCAAGGATTCTTACTTTGCCTTCATGGCCAGTGCTATGCACTGGGTGAAGCAAAGACCTGGACATGGCCTGGAATGGATAGGATCTTTTACTATGTACAGTGATGCTACTGAGTACAGTGAAAACTTCAAGGGCAAGGCCACATTGACTGCAAACACATCCTCGAGCACAGCCTACATGGAACTCAGCAGCCTCACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1015560,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.50.143*01",
        "confirmed": True,
        "end": 951330,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.50.143",
        "id": "J558.50.143*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-50*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-50",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCAGTGGGTAAAACAGAGGCCTGGACAGGGCCTTGAGTGGATCGGAGAGATTGATCCTTCTGATAGCTATACTAACTACAATCAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 951036,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.51pg.144*01",
        "confirmed": True,
        "end": 939709,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.51pg.144",
        "id": "J558.51pg.144*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-51*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-51",
        "length": 300,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTTCCTCAGTCTGGTTCTGAGGTGGGGCGGACTGGTGCCTCAGTGAAGATGTTCTGCAAGGCTCCTGGCTACACATTCACTAACTACTATATGTATTGATTAAAGCAGAGTCATGGAGATAGCCTAGAGTGGATTTGATATATTTATCCTGGAAATGGTCTTACTAGCTATGCCAAGAAGTTCAAAGGCAAGGCCACATTGACTATAGACAATTCAGCCAGCACAGCCTACATGCAGCTCAGCAGCATGACATCTGAAGCCTCTGATGACTATTGTTGTGCTAGACAAGTG",
            "read frame": 0,
            "strand": True
        },
        "start": 939409,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.52.145*01",
        "confirmed": True,
        "end": 925591,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.52.145",
        "id": "J558.52.145*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-52*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-52",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGTCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCATTGGGTGAAGCAGAGGCCTATACAAGGCCTTGAATGGATTGGTAACATTGACCCTTCTGATAGTGAAACTCACTACAATCAAAAGTTCAAGGACAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 925297,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.53.146*01",
        "confirmed": True,
        "end": 912675,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.53.146",
        "id": "J558.53.146*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-53*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-53",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAAATATTAATCCTAGCAATGGTGGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 912381,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.54.148*01",
        "confirmed": True,
        "end": 877403,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.54.148",
        "id": "J558.54.148*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-54*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-54",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTAAGGCCTGGGACTTCAGTGAAGGTGTCCTGCAAGGCTTCTGGATACGCCTTCACTAATTACTTGATAGAGTGGGTAAAGCAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAGTGATTAATCCTGGAAGTGGTGGTACTAACTACAATGAGAAGTTCAAGGGCAAGGCAACACTGACTGCAGACAAATCCTCCAGCACTGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 877109,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.55.149*01",
        "confirmed": True,
        "end": 862983,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.55.149",
        "id": "J558.55.149*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-55*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-55",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATAACCTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAGATATTTATCCTGGTAGTGGTAGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 862689,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.56.150*01",
        "confirmed": True,
        "end": 828276,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.56.150",
        "id": "J558.56.150*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-56*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-56",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTCCTGGCTATACCTTCACCAGCCACTGGATGCAGTGGGTAAGACAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAGAGATTTTTCCTGGAAGTGGTAGTACTTATTATAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 827982,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.57pg.152*01",
        "confirmed": True,
        "end": 802654,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.57pg.152",
        "id": "J558.57pg.152*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-57*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-57",
        "length": 292,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCTAGCTGCATCAATCTGGGGCTGAGCTGGTGAAGCCTGGGATCTCTGTGAAGATGTTATACATGGTTTCTGGTTATACATTTACTGAATACTACATGCCCTGGTCAAGCAGAATCATGGAAAGATCCTTGAGTGGATTGGAAATGTTAATATTTACAATTGTGGTATAACTACAATAAAAATTTCAAGGACAAGGACACATCAACTGTAGACTATTCCTCCAGTACAGCCTATATGTTGCTTGGCAAAGTGACATCTGAGGATTCTAAGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 802362,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.58.154*01",
        "confirmed": True,
        "end": 758912,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.58.154",
        "id": "J558.58.154*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-58*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-58",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTTCAGCAGTCTGGAGCTGAGCTGGTGAGGCCTGGGTCCTCAGTGAAGATGTCCTGCAAGACTTCTGGATATACATTCACAAGCTACGGTATAAACTGGGTGAAGCAGAGGCCTGGACAGGGCCTGGAATGGATTGGATATATTTATATTGGAAATGGTTATACTGAGTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTTCAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCAATCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 758618,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.59.155*01",
        "confirmed": True,
        "end": 735996,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.59.155",
        "id": "J558.59.155*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-59*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-59",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGACTTCAGTGAAGTTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTAAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATCGGAGTGATTGATCCTTCTGATAGTTATACTAACTACAATCAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 735702,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.5pg.94*01",
        "confirmed": True,
        "end": 1515499,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.5pg.94",
        "id": "J558.5pg.94*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-8*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-8",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGATCAGCTGCATCAGTCTGAAGCTGAGCTGCAGCAACCTGGGACATCAGTGAAGATGCCCTGAAAGGCTACTGGCTACACCTTCACTAAGTATCGAATGTGTTGGGTGAGGCAGAAGCTTGGACAGGGCCTGGAATGGATTGCATCTGTTGATCCTGGAAATAGTAATACTGAATACAATCAGAAGTTCAAAGGCAAGGCCACACTAACAGAACACAAATCCTCCAGCACAGCCTACATAGAGCTTAGCAACCTGACCTCTGAGGACTCTGCTGTCTATTACTGTACAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1515205,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.6.96*01",
        "confirmed": True,
        "end": 1487539,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.6.96",
        "id": "J558.6.96*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-9*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-9",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAGCTGCAGCAGTCTGGAGCTGAGCTGATGAAGCCTGGGGCCTCAGTGAAGCTTTCCTGCAAGGCTACTGGCTACACATTCACTGGCTACTGGATAGAGTGGGTAAAGCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGAGATTTTACCTGGAAGTGGTAGTACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACATTCACTGCAGATACATCCTCCAACACAGCCTACATGCAACTCAGCAGCCTGACAACTGAGGACTCTGCCATCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1487245,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.60pg.156*01",
        "confirmed": True,
        "end": 726547,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.60pg.156",
        "id": "J558.60pg.156*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-60*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-60",
        "length": 300,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTTCCTCAGTCTGGATCTGAGGTGGGGCGGACTGGTGCCTCAGTGAAGATGTTCTGCAAGGCTTCTGGCTACACATTCACTAACTACTATATGCATTGATTAAAGCAGAGTCATGGAGATAGCCTAGAGTGGATTTGATATATTTATCCTGGAAATGGTCTTACTAGCTATGCCAAGAAGTTCAAAGGCAAGGCCACATTGACTATAGACAATTCAGCCAGCACAGCCTACATGCAGCTCAGCAGCATGACATCTGAAGCCTCTGATGACTATTACTGTGCTAGACAAGTG",
            "read frame": 0,
            "strand": True
        },
        "start": 726247,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.61.157*01",
        "confirmed": True,
        "end": 711938,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.61.157",
        "id": "J558.61.157*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-61*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-61",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGTCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGGATTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAATGGATTGGTAACATTTACCCTTCTGATAGTGAAACTCACTACAATCAAAAGTTCAAGGACAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 711644,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.62pg.158*01",
        "confirmed": True,
        "end": 699092,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.62pg.158",
        "id": "J558.62pg.158*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-62*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-62",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCGCAAGGCTTCTGGCTACACTTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGAGAAGGCCTTGAGTGGATTGGAAATATTTATCCTGGTAGTAGTAGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 698799,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.63pg.159*01",
        "confirmed": True,
        "end": 684383,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.63pg.159",
        "id": "J558.63pg.159*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-62-1*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-62-1",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCAGTGGGTAAAACAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAGAGATTTTTCCTGGAAGTGGTAGTACTTATTATAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACAGCAGAGAACTCTGCAATCTATCTGTGCAAGGAA",
            "read frame": 0,
            "strand": True
        },
        "start": 684089,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.64.162*01",
        "confirmed": True,
        "end": 624662,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.64.162",
        "id": "J558.64.162*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-62-2*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-62-2",
        "length": 302,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTGAAACCCGGGGCATCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACTGAGTATACTATACACTGGGTAAAGCAGAGGTCTGGACAGGGTCTTGAGTGGATTGGGTGGTTTTACCCTGGAAGTGGTAGTATAAAGTACAATGAGAAATTCAAGGACAAGGCCACATTGACTGCGGACAAATCCTCCAGCACAGTCTATATGGAGCTTAGTAGATTGACATCTGAAGACTCTGCGGTCTATTTCTGTGCAAGACACGAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 624360,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.65.163*01",
        "region": "VH",
        "end": 610072,
        "family": "J558",
        "framed": True,
        "functionality": "O",
        "gene": "J558.65.163",
        "id": "J558.65.163*01-C57BL/6",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-62-3",
        "length": 294,
        "read frame": 0,
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCAGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACGAGGCCTTGAGTGGATTGGAAGGATTGATCCTAATAGTGGTGGTACTAAGTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAACCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 609778,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.66.165*01",
        "confirmed": True,
        "end": 575445,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.66.165",
        "id": "J558.66.165*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-63*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-63",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTAAGGCCTGGGACTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACCTTCACTAACTACTGGATAGGTTGGGCAAAGCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGATATTTACCCTGGAGGTGGTTATACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGTTCAGCAGCCTGACATCTGAGGACTCTGCCATCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 575151,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.67.166*01",
        "confirmed": True,
        "end": 563533,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.67.166",
        "id": "J558.67.166*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-64*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-64",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGTTGTCCTGCAAGGCTTCTGGCTACACTTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAATGATTCATCCTAATAGTGGTAGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAACTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 563239,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.68pg.168*01",
        "confirmed": True,
        "end": 539022,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.68pg.168",
        "id": "J558.68pg.168*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-65*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-65",
        "length": 292,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCTAGCTGCAGCAGTCAGGGGCTGAGCTGGTGAAGCCTGATGGCTCTGTGAAGATGTCATGCAAGGCTTCTGGTTATACATTTACTGAATACCACATGCTCTAGTCAAGCAGAATCATGGAAAGACCCTTGAATGGATTTTCAATATTAATACTTAAAATGGTGGTATAACTACAGTGAAAATTTCAAGGGCAAGGGTACATTAACTGTAGACAAATCCTCCAGCACACCCTATTTGTTGCTTAGCAAATTGACATCTGAGGATTCTGTGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 538730,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.69.170*01",
        "confirmed": True,
        "end": 477968,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.69.170",
        "id": "J558.69.170*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-66*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-66",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACAGCTTCACAAGCTACTATATACACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTATCCTGGAAGTGGTAATACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACGGCAGACACATCCTCCAGCACTGCCTACATGCAGCTCAGCAGCCTAACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 477674,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.70pg.171*01",
        "confirmed": True,
        "end": 467138,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.70pg.171",
        "id": "J558.70pg.171*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-67*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-67",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGGCCTGAGCTGGTGAGGCCTGGGGTCTCAGTGAAGATTTCCTGCAAGGGTTCCGGCTACACATTCACTGATTATGCTATGCACTGGGTGAAACAGAGTCATGCAAAGAGTCTAGAGTGGATTGGAGTTATTAGTACTTACTATGGTGATGCTAGCTACAACCAGAAGTTCAAGGACAAGGCCACAATGACTGTAGACAAATCCTCCAGCACAGCCTATATGGAACTTGCCAGACTGACATCTGAGGACTCTGCCGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 466844,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.71pg.172*01",
        "confirmed": True,
        "end": 459738,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.71pg.172",
        "id": "J558.71pg.172*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-68*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-68",
        "length": 284,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CTGGTTCAGCTGCAGCAATGTGAAGCTGAGGTTGTGAGACTCAGTGATGGTGTCCTGCTCAGCTTCTGTCTACATATTTAGCAACTACTATGAAGTGTATAAAGCAGAGGCCTGGACAGGCTCTTGAGTGGATTCGATGAATTGATTCTTGAAGTTTTGGTACTACCTATAATCAGAAGTTCAAAGGCCAGGCCACATCTACTGTAGACAAATCCTCCATCACAGCCTACCTGCAACTAATATCCTGACATCTGAGAACTCTGAAGACTATTACTGTGAAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 459454,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.72.173*01",
        "confirmed": True,
        "end": 447917,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.72.173",
        "id": "J558.72.173*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-69*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-69",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGATGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATCGGAGAGATTGATCCTTCTGATAGTTATACTAACTACAATCAAAAGTTCAAGGGCAAGTCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 447623,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.73pg.175*01",
        "confirmed": True,
        "end": 379189,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.73pg.175",
        "id": "J558.73pg.175*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-70*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-70",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAGGCCTGGGACTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTATACCTTCTTCACCTACTGGATGAACTGGGTGTAGCAGAGGCCTGGACAGGGCCTTGAGTGGATTGGACAGATTTTTCCTGCAAGTGGTAGTACTTACTACAATGAGATGTACAAGGACAAGGCCGCATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAAGACACTGCTGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 378895,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.74.176*01",
        "confirmed": True,
        "end": 328873,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.74.176",
        "id": "J558.74.176*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-71*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-71",
        "length": 302,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTGAAACCCGGGGCATCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACTGAGTATACTATACACTGGGTAAAGCAGAGGTCTGGACAGGGTCTTGAGTGGATTGGGTGGTTTTACCCTGGAAGTGGTAGTATAAAGTACAATGAGAAATTCAAGGACAAGGCCACATTGACTGCGGACAAATCCTCCAGCACAGTCTATATGGAGCTTAGTAGATTGACATCTGAAGACTCTGCGGTCTATTTCTGTGCAAGACACGAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 328571,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.75.177*01",
        "confirmed": True,
        "end": 313094,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.75.177",
        "id": "J558.75.177*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-72*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-72",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACGAGGCCTTGAGTGGATTGGAAGGATTGATCCTAATAGTGGTGGTACTAAGTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAACCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 312800,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.76pg.179*01",
        "confirmed": True,
        "end": 281464,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.76pg.179",
        "id": "J558.76pg.179*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-73*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-73",
        "length": 295,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTAAGGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTCACTAACTACTGGATAGGTTGGGTAAAGCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGATATTTACCCTGGAGATGGTTATACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGTGCATGCAGCCTGACCTCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 281169,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.77.180*01",
        "confirmed": True,
        "end": 269407,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.77.180",
        "id": "J558.77.180*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-74*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-74",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAACTGCAGCAGCCTGGGGCTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGGTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGCCAAGGCCTTGAGTGGATTGGAAGGATTCATCCTTCTGATAGTGATACTAACTACAATCAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAATA",
            "read frame": 0,
            "strand": True
        },
        "start": 269113,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.78.182*01",
        "confirmed": True,
        "end": 238104,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.78.182",
        "id": "J558.78.182*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-75*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-75",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTACAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACACCTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTTTCCTGGAAGTGGTAGTACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTTACTGTAGACAAATCCTCCAGCACAGCCTACATGTTGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 237810,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.79.184*01",
        "confirmed": True,
        "end": 224171,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.79.184",
        "id": "J558.79.184*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-76*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-76",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGAAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACTTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGCAAGGATTTATCCTGGAAGTGGTAATACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGAAAAATCCTCCAGCACTGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCTGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 223877,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.7pg.97*01",
        "confirmed": True,
        "end": 1473229,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.7pg.97",
        "id": "J558.7pg.97*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-10*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-10",
        "length": 308,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CCAGCCCAGCTGCAGCATTCTGGGGATGACATGGAGGAGCCTGGGTCCTCAGTGAAGTTTTCCTGCATGGCTTCTGGCTATACCTTCACTGACCACTCTATGCATGTGGTAAAACAATGAAACAGAGGCCTGGACAGGGACTGGAGTGGATTGAATGGGTTGGACCTACATGTGGTGGTACTGTATATGCTAGGAAGTTTCAAGGCAAGGCCACACTGACTGTAGACAAATCAGCCATCACAACCTACATGCAGCTCAGTAGCCTGACATCTGAAAACTCTGCAGTCTATTACTGTGCCATGCAAGGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1472921,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.8.98*01",
        "confirmed": True,
        "end": 1458865,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.8.98",
        "id": "J558.8.98*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-11*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-11",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATCCAGCTGCAACAGTCAGGAGCTGAGCTGGCGAGTCCTGGGGCATCAGTGACACTGTCCTGCAAGGCTTCTGGCTACACATTTACTGACCATATTATGAATTGGGTAAAAAAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAAGGATTTATCCAGTAAGTGGTGAAACTAACTACAATCAAAAGTTCATGGGCAAGGCCACATTCTCTGTAGACCGGTCCTCCAGCACAGTGTACATGGTGTTGAACAGCCTGACATCTGAGGACCCTGCTGTCTATTACTGTGGAAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1458571,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.80.186*01",
        "confirmed": True,
        "end": 210185,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.80.186",
        "id": "J558.80.186*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-77*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-77",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGAAGCAGTCTGGAGCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACACCTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAAAGATTGGTCCTGGAAGTGGTAGTACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 209891,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.81.187*01",
        "confirmed": True,
        "end": 203279,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.81.187",
        "id": "J558.81.187*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-78*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-78",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAGCTGCAACAGTCTGACGCTGAGTTGGTGAAACCTGGAGCTTCAGTGAAGATATCCTGCAAGGTTTCTGGCTACACCTTCACTGACCATACTATTCACTGGATGAAGCAGAGGCCTGAACAGGGCCTGGAATGGATTGGATATATTTATCCTAGAGATGGTAGTACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAACAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 202985,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.82pg.188*01",
        "confirmed": True,
        "end": 183956,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.82pg.188",
        "id": "J558.82pg.188*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-79*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-79",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAGCTGCAGTAGTCTGGAGCTGAACTGGTGAAACCAGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACAAACTACATCATAAACTGGGTAAAGGAAGGGCCTGGACATGGCCTTGAGTGGATTGGATGGATTTCTCCTGAATATGGTCATACTTACTACAATCAGAAGTTCAAGGGCAAGGCCACATTCACTGCAGACACATCCTCCAGCACAGCCTACATGGAGCTCAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 183662,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.83.189*01",
        "confirmed": True,
        "end": 159953,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.83.189",
        "id": "J558.83.189*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-80*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-80",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAGCTGCAGCAGTCTGGGGCTGAGCTGGTGAAGCCTGGGGCCTCAGTGAAGATTTCCTGCAAAGCTTCTGGCTACGCATTCAGTAGCTACTGGATGAACTGGGTGAAGCAGAGGCCTGGAAAGGGTCTTGAGTGGATTGGACAGATTTATCCTGGAGATGGTGATACTAACTACAACGGAAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 159659,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.84.190*01",
        "confirmed": True,
        "end": 152016,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.84.190",
        "id": "J558.84.190*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-81*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-81",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAGCTGCAGCAGTCTGGAGCTGAGCTGGCGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACAAGCTATGGTATAAGCTGGGTGAAGCAGAGAACTGGACAGGGCCTTGAGTGGATTGGAGAGATTTATCCTAGAAGTGGTAATACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCGTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 151722,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.85.191*01",
        "confirmed": True,
        "end": 119766,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.85.191",
        "id": "J558.85.191*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-82*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-82",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCCTCAGTGAAGATTTCCTGCAAGGCTTCTGGCTACGCATTCAGTAGCTCCTGGATGAACTGGGTGAAGCAGAGGCCTGGAAAGGGTCTTGAGTGGATTGGACGGATTTATCCTGGAGATGGAGATACTAACTACAATGGGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAACTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTACTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 119472,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.86.192*01",
        "confirmed": True,
        "end": 108481,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.86.192",
        "id": "J558.86.192*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-83*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-83",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACATTCACTGACTACTACATGCACTGGGTGAAGCAGAAGCCTGGGAAGGGCCTTGAGTGGATTGGAGAGATTTATCCTGGAAGTGGTAATACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 108187,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.87.193*01",
        "confirmed": True,
        "end": 91596,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.87.193",
        "id": "J558.87.193*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-84*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-84",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACACCTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTATCCTGGAAGCGGTAATACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 91302,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.88.194*01",
        "confirmed": True,
        "end": 72270,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.88.194",
        "id": "J558.88.194*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-85*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-85",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGTTGTCCTGCAAGGCTTCTGGCTACACCTTCACAAGCTACGATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTATCCTAGAGATGGTAGTACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCGTACATGGAGCTCCACAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 71976,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.89pg.195*01",
        "confirmed": True,
        "end": 62639,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.89pg.195",
        "id": "J558.89pg.195*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-86*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-86",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCCTCAGTGAAGATTTCCTGCAAGGCTTTTGGCTACACCTTCACAAACCATCATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGCCTGGACTGGATTGGATATATTAATCCTTATAATGATTATACTAGCTACAGAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTATATGGAGCTTAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 62343,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.9.99*01",
        "confirmed": True,
        "end": 1455258,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.9.99",
        "id": "J558.9.99*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-12*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-12",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGCTTATCTACAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCCTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACATTTACCAGTTACAATATGCACTGGGTAAAGCAGACACCTAGACAGGGCCTGGAATGGATTGGAGCTATTTATCCAGGAAATGGTGATACTTCCTACAATCAGAAGTTCAAGGGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAAGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1454964,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.16.76*01",
        "confirmed": True,
        "end": 1719404,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "PG.16.76",
        "id": "PG.16.76*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-1*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-1",
        "length": 291,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAAATCCAACTGCAACCGTCTGGGGCCGACCTGATTAATCCTGGATCCTCAATGAAGGTGTCTTGCAAGGTTCCTGGCTACAGTTTCACTAGCTACTATATGACCTGGGTAAAAGAGAGTCCTGGACAGGGTCTAGAATAGATTAGAGAAACCACCCTAGTAGTTGCAGTATAAGTTATGCACAGAAGTTTCAAGGACACATCTCCATGACTAGGCACATATCCTCCAGCATGGCCTACATGGAGCTCAGCAGGCTGACCTCTGAGGACACTTCTGTTTATTGCTGTTTAT",
            "read frame": 0,
            "strand": True
        },
        "start": 1719113,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.17.87*01",
        "confirmed": True,
        "end": 1589494,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "PG.17.87",
        "id": "PG.17.87*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-3*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-3",
        "length": 302,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "AAGGTCCAGCTGCAGCAGTCTGGAACTGAGTTTGGGAGGTCTGGAGCCTCAGAGAAGTAGTCCTGCAAGTCTTCTGGTTACAAATTCACCAACTATTATATGCACCGGATAGGGCAGAATCGTGGAAATAGCCTAGACTAGATTGGATATATTTATCGTGGAATTGGCCATCCTGGCTATGCCCACAACTTCAGGGATAAGACCACACTGACTTCAGACAAATCCTCTAGCACAGCCTACATGGAATTTCACATCCATTACTCCAGGAATCTTCTTGGGGTATCTCCCCATGACGTAACTCA",
            "read frame": 0,
            "strand": True
        },
        "start": 1589192,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.18.92*01",
        "confirmed": True,
        "end": 1538136,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "PG.18.92",
        "id": "PG.18.92*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-6*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-6",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGAACTGAGTTTGGGTGGACTAGAGCCTCAGGGAAGTAGTCCTGCAAGTCTTCTGGTTACAAATTCACCAACTACTATATGCACTGGATAGGGCAGAATCGTGGAAATAGCCTAGACTAGATTGGATATATTTATCGTGGAATTGGCCATCCTGGCTATGCTCACAACTTCAAGGACAAGACCACACTGACTTCAGACAAATCCTCTAGCACAGCCTACATGGAATTTCACATCCATTACTCCAGGAATCTTCTTGGGGTATCTCCCCATGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1537843,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "7183.5pg.7*01",
        "confirmed": True,
        "end": 2467796,
        "family": "J606",
        "framed": True,
        "functionality": "P",
        "gene": "7183.5pg.7",
        "id": "7183.5pg.7*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-1*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-1",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTAGTGGAGTCTAGGGGAGGCTTAATGCAGCCTGAAAGGTCCGTGATACTCTCCTGTGCTGCCTCTGGATGCACTGTCAGTGACTACTGGTTGGCCTGGGTTTGCCAGGCTCCAAAGAAGGGGCTGGATGGGGGGGGGCATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGAAGAAGTGGGCTGGACATAACATATTCAGAATCAATATTTAGAGATTCTAATCCTTAAAGACATCACTTTGACCAAGTATATATGCAACATGTTATGGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 2467502,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J606.1.79*01",
        "confirmed": True,
        "end": 1679397,
        "family": "J606",
        "framed": True,
        "functionality": "F",
        "gene": "J606.1.79",
        "id": "J606.1.79*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-3*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-3",
        "length": 299,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGAAGCTTGAGGAGTCTGGAGGAGGCTTGGTGCAACCTGGAGGATCCATGAAACTCTCCTGTGTTGCCTCTGGATTCACTTTCAGTAACTACTGGATGAACTGGGTCCGCCAGTCTCCAGAGAAGGGGCTTGAGTGGGTTGCTCAAATTAGATTGAAATCTGATAATTATGCAACACATTATGCGGAGTCTGTGAAAGGGAGGTTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTCTACCTGCAAATGAACAACTTAAGGGCTGAAGACACTGGAATTTATTACTGCACAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1679098,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J606.2.80*01",
        "confirmed": True,
        "end": 1664635,
        "family": "J606",
        "framed": True,
        "functionality": "F",
        "gene": "J606.2.80",
        "id": "J606.2.80*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-4*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-4",
        "length": 300,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GATGTGAACCTGGAAGTGTCTGGAGGAGGCTTAGTTAAACCTGGAGGATCCATGCAACTCTTTTGTGTAGCCTCTGGATTTACTTTTGTAGATGGCTGGATGGACTGGGTCCGCCAGTCTCCAGAGAAGGGTCTTGAGTGGGTTGCTGAAATTGCAAACAAAGCTAATAATTATGCAACATATTATCCCGAGTCTGTGAAAGGCAGATTCACCATCTCAAGAGATGATTTCAAAAGTAGTGTCTACCTGCATATGAACAGCTTAAGAGCCGAAGATACAGGCATTTATTACTGTACAAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1664335,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J606.3.81*01",
        "confirmed": True,
        "end": 1654513,
        "family": "J606",
        "framed": True,
        "functionality": "F",
        "gene": "J606.3.81",
        "id": "J606.3.81*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-5*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-5",
        "length": 300,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGAAAATTGAGGAGTCAGGAGGAGGCTTGGTCCAACCTGGAGGATCCATGAAACTCTCTTGTGCAGCCTCTGGATTCACTTTCAGTGATTACAGGATGGACTGGGTCCACCACTCTACAGAGAATGGGTTGGAGTGGGTTGCTGAAATTAGAAACAAAGCTAGTAATTATGCAACATATTATGTGGAGTCTGTGAATGGGAGGTTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTCTACCTGCAAATGAACAGCTTAAGAGCTGAAGATACTGGCATTTATTACTGTACAAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1654213,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J606.4.82*01",
        "confirmed": True,
        "end": 1636319,
        "family": "J606",
        "framed": True,
        "functionality": "F",
        "gene": "J606.4.82",
        "id": "J606.4.82*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-6*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-6",
        "length": 300,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGAAGCTTGAGGAGTCTGGAGGAGGCTTGGTGCAACCTGGAGGATCCATGAAACTCTCTTGTGCTGCCTCTGGATTCACTTTTAGTGACGCCTGGATGGACTGGGTCCGCCAGTCTCCAGAGAAGGGGCTTGAGTGGGTTGCTGAAATTAGAAACAAAGCTAATAATCATGCAACATACTATGCTGAGTCTGTGAAAGGGAGGTTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTCTACCTGCAAATGAACAGCTTAAGAGCTGAAGACACTGGCATTTATTACTGTACCAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1636019,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "J606.5.83*01",
        "confirmed": True,
        "end": 1615483,
        "family": "J606",
        "framed": True,
        "functionality": "F",
        "gene": "J606.5.83",
        "id": "J606.5.83*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-7*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-7",
        "length": 300,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGAGAAGCTGGATGAGTCTGGAGGAGGCTTGGTGCAACCTGGGAGGTCCATGAAACTCTCCTGTGTTGCCTCTGGATTCACTTTTACTAACTCCTGGATGAACTGGTTCTGCCAGTCTCCAGAGAAAGGACTGGAGTGGGTAGCACAAATTAAAAGCAAACCTTATAATTATGAAACATATTATTCAGATTCTGTGAAAGGCAGATTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTATACCTGCAAATGAACAACTTAAGAGCTGAAGACACGGGCATCTATTACTGTACATGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1615183,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.10.56*01",
        "confirmed": True,
        "end": 1981725,
        "family": "J606",
        "framed": True,
        "functionality": "P",
        "gene": "PG.10.56",
        "id": "PG.10.56*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-2*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-2",
        "length": 295,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGTTGCAGCTTATGGAGTCTTGGGGAAGCTTGTTAAAGCTCCAGGGTTCTGTAAGACATTCTTGTGCAGCCGCTGGATTCACTTTCAGAGAATACTATATCAGATGGGTCCTGTAGATTTTAGAAAAAGATCTTGAGTCCTTGATTTAATTACAAACACAGCTGGGGGGTGATTACAGAGTATGCTTCCTATGTGAGAAGGGCACTTCCCATTTCAAGAGATCATACAAAAAATAAAATTTCATGCCTCCAAATAAACAGATTGGAGACTTCTTTTTTCCCTGAGAAGAAATGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1981430,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.10.33*01",
        "confirmed": True,
        "end": 2275011,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.10.33",
        "id": "Q52.10.33*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-6-8*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-6-8",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACTGTCTCTGGGTTCTCATTAACCAGCTATGGTGTAGACTGGGTTCGCCAGTCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTAATATGGGGTGTTGGAAGCACAAATTATAATTCAGCTCTCAAATCCAGACTGAGCATCAGCAAGGACAACTCCAAGAGCCAAGTTTTCTTAAAAATGAACAGTCTGCAAACTGATGACACAGCCATGTACTACTGTGCCAGTGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2274718,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.11.34*01",
        "confirmed": True,
        "end": 2263835,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.11.34",
        "id": "Q52.11.34*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-7*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-7",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAAGTGCAGATGAAGGAGTCAGGACCTGACCTTGTGCAGCCATCACAGACTCTGTCTCTCACCTGCACTGTCTCTGGGTTCTCATTAAGTAGCTATGGTGTACATTGGTTTCGCAAGCCTCCGAGAAAGGGATTGGAATGGTTGGGAGGAATATGGTCTGGTGGAAGCATATACTATACTCCAGCTCTCAGTTCCCGACTGAGTGTCAGCAGGGACACCTCTAAGAGCCAAGTTTTCTTTAAAATGAGCAGTCTGCAAAGTGAAGACACGGCTGTGTACCACTGTGCCAGATA",
            "read frame": 0,
            "strand": True
        },
        "start": 2263542,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.12pg.39*01",
        "confirmed": True,
        "end": 2193788,
        "family": "Q52",
        "framed": True,
        "functionality": "P",
        "gene": "Q52.12pg.39",
        "id": "Q52.12pg.39*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-8*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-8",
        "length": 299,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ACCCTGTAGAGGTACAACTGAAGTAGTCATGACCTGACCTGCTGCAATCCTCACAGGCTGTCTCTCATCTGCACTGTCTCTGGCTTCTCATTAACCATCTATGGTGTAAACTGGGTCCTCCAGCCACTAGGAAGGGGATTGCAGAGGATGCCAGCAATATGGAGAGGTGAAAGCACAGAGTATAATTCAGATCTCAAATCCTGAATCAGCATCAGTAGGGACACATCCAAGAGTCAACTTTTCTTAAAACTGAACAGAATTCAAACTGAGGACATAGCCATCTATGACCCTGTCAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2193489,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.13.40*01",
        "confirmed": True,
        "end": 2192053,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.13.40",
        "id": "Q52.13.40*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-9*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-9",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACTTGCACTGTCTCTGGGTTTTCATTAACCAGCTATGGTGTAGACTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTAATATGGGGTGGTGGAAGCACAAATTATAATTCAGCTCTCATGTCCAGACTGAGCATCAGCAAAGACAACTCCAAGAGCCAAGTTTTCTTAAAAATGAACAGTCTGCAAACTGATGACACAGCCATGTACTACTGTGCCAAACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2191760,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.1pg.2*01",
        "confirmed": True,
        "end": 2497009,
        "family": "Q52",
        "framed": True,
        "functionality": "P",
        "gene": "Q52.1pg.2",
        "id": "Q52.1pg.2*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-1*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-1",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGTTCCAGCTGAAGCAGTCAGGACCTGCCCTTGTGCCGCCTTCACGGAGCCTGTCCATCACATGCACTGTCTCTGGGTTCTCATTAACCAGCTGTGGTATAAGCTGGCTTCGATGGTCTCCAGGAAACGGAGATAGTGGAATTTTACAACTGTCCCATAGGCCTAGAAAACCCAGTCATCACCCCTGTTCCTCTTCTGCATCCTTACCAATACGTGATCCCCAATCCTTACCCACTATTCATCAAGTATATTATATTTTGCTACTGTTCCCCATTTATACAGGACTCGCCAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2496715,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.2.4*01",
        "confirmed": True,
        "end": 2482989,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.2.4",
        "id": "Q52.2.4*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-2*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-2",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTGCAGCTGAAGCAGTCAGGACCTGGCCTAGTGCAGCCCTCACAGAGCCTGTCCATCACCTGCACAGTCTCTGGTTTCTCATTAACTAGCTATGGTGTACACTGGGTTCGCCAGTCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTGATATGGAGTGGTGGAAGCACAGACTATAATGCAGCTTTCATATCCAGACTGAGCATCAGCAAGGACAATTCCAAGAGCCAAGTTTTCTTTAAAATGAACAGTCTGCAAGCTGATGACACAGCCATATATTACTGTGCCAGAAA",
            "read frame": 0,
            "strand": True
        },
        "start": 2482696,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.3.8*01",
        "confirmed": True,
        "end": 2460072,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.3.8",
        "id": "Q52.3.8*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-3*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-3",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACTGTCTCAGGGTTCTCATTAACCAGCTATGGTGTAAGCTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTAATATGGGGTGACGGGAGCACAAATTATCATTCAGCTCTCATATCCAGACTGAGCATCAGCAAGGATAACTCCAAGAGCCAAGTTTTCTTAAAACTGAACAGTCTGCAAACTGATGACACAGCCACGTACTACTGTGCCAAACC",
            "read frame": 0,
            "strand": True
        },
        "start": 2459779,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.4pg.12*01",
        "end": 2423778,
        "family": "Q52",
        "framed": True,
        "functionality": "P",
        "gene": "Q52.4pg.12",
        "id": "Q52.4pg.12*01-C57BL/6",
        "length": 288,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ACTTGCATATAAAAGAGTGAGGACTTGAACTGGTGTAGCCTTCACAGACCCTGCCCCATACCTGCACTGTCTCTGGCTTCTCATTATTCCAGCTACCATGTGCACTGAGTGTGTCAGTCTACAGCAAAGGGTCCACAGTGGATGGGAGCAATGTGTAGTGGGGAAACACAGCATACGGTTCAGCTCTCAACTCCCAACTCGTCATCAATAGGGACACATCTATGGGCAAGTGTTCTTAAAACTGTACAATCTTCAAAAAGGGAAACAGTGATGTACTACTGTGCCAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2423490,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.5.13*01",
        "confirmed": True,
        "end": 2417965,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.5.13",
        "id": "Q52.5.13*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-4*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-4",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTGCAGCTGAAGCAGTCAGGACCTGGCCTAGTGCAGCCCTCACAGAGCCTGTCCATCACCTGCACAGTCTCTGGTTTCTCATTAACTAGCTATGGTGTACACTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTGATATGGAGTGGTGGAAGCACAGACTATAATGCTGCTTTCATATCCAGACTGAGCATCAGCAAGGACAACTCCAAGAGCCAAGTTTTCTTTAAAATGAACAGTCTGCAAGCTGATGACACTGCCATATACTACTGTGCCAAAAA",
            "read frame": 0,
            "strand": True
        },
        "start": 2417672,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.6pg.17*01",
        "end": 2388792,
        "family": "Q52",
        "framed": True,
        "functionality": "P",
        "gene": "Q52.6pg.17",
        "id": "Q52.6pg.17*01-C57BL/6",
        "length": 300,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATCGTCTCCTACTTGCATCTGAAAGAGCGAGGACTTGAACTGGTGTAGCCTTCACAGACCCTGCCCCATACCTGCACTGTCTCTGGCTTCTCATTATTCAAGCTACCATGTGCACTGAGTGTGTCAGTCTACAGCAAAGGGTCCACAGTGGATGGGAGCAATGTGTAGTGGGGAAACACAGCATACGGTTCAGCTCTCAACTCCCAACTCGTCATCAATAGGGACACATCTATGGGCAAGTGTTCTTAAAACTGAACAATCTTCAAAAAGGAAAACAGTGATGTACTACTGTGCCAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2388492,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.7.18*01",
        "confirmed": True,
        "end": 2385774,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.7.18",
        "id": "Q52.7.18*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-5*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-5",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTGCAGCTGAAGCAGTCAGGACCTGGCCTAGTGCAGCCCTCACAGAGCCTGTCCATAACCTGCACAGTCTCTGGTTTCTCATTAACTAGCTATGGTGTACACTGGGTTCGCCAGTCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTGATATGGAGAGGTGGAAGCACAGACTACAATGCAGCTTTCATGTCCAGACTGAGCATCACCAAGGACAACTCCAAGAGCCAAGTTTTCTTTAAAATGAACAGTCTGCAAGCTGATGACACTGCCATATACTACTGTGCCAAAAA",
            "read frame": 0,
            "strand": True
        },
        "start": 2385481,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.8.22*01",
        "confirmed": True,
        "end": 2354482,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.8.22",
        "id": "Q52.8.22*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-6*03",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-6",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACCGTCTCAGGGTTCTCATTAACCAGCTATGGTGTACACTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGTAGTGATATGGAGTGATGGAAGCACAACCTATAATTCAGCTCTCAAATCCAGACTGAGCATCAGCAAGGACAACTCCAAGAGCCAAGTTTTCTTAAAAATGAACAGTCTCCAAACTGATGACACAGCCATGTACTACTGTGCCAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2354189,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "Q52.9.29*01",
        "confirmed": True,
        "end": 2301298,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.9.29",
        "id": "Q52.9.29*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-9-1*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-9-1",
        "length": 293,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACTGTCTCTGGGTTCTCATTAACCAGCTATGCTATAAGCTGGGTTCGCCAGCCACCAGGAAAGGGTCTGGAGTGGCTTGGAGTAATATGGACTGGTGGAGGCACAAATTATAATTCAGCTCTCAAATCCAGACTGAGCATCAGCAAAGACAACTCCAAGAGTCAAGTTTTCTTAAAAATGAACAGTCTGCAAACTGATGACACAGCCAGGTACTACTGTGCCAGAAA",
            "read frame": 0,
            "strand": True
        },
        "start": 2301005,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "S107.1.42*01",
        "confirmed": True,
        "end": 2174739,
        "family": "S107",
        "framed": True,
        "functionality": "P",
        "gene": "S107.1.42",
        "id": "S107.1.42*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV7-1*01",
        "imgt family name": "IGHV7",
        "imgt gene name": "IGHV7-1",
        "length": 304,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTGGAATCTGGAGGAGGCTTGGTACAGTCTGGGCGTTCTCTGAGACTCTCCTGTGCAACTTCTGGGTTCACCTTCAGTGATTTCTACATGGAGTGGGTCCGCCAAGCTCCAGGGAAGGGACTGGAGTGGATTGCTGCAAGTAGAAACAAAGCTAATGATTATACAACAGAGTACAGTGCATCTGTGAAGGGTCGGTTCATCGTCTCCAGAGACACTTCCCAAAGCATCCTCTACCTTCAGATGAATGCCCTGAGAGCTGAGGACACTGCCATTTATTACTGTGCAAGAGATG",
            "read frame": 0,
            "strand": True
        },
        "start": 2174435,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "S107.2pg.43*01",
        "confirmed": True,
        "end": 2159124,
        "family": "S107",
        "framed": True,
        "functionality": "F",
        "gene": "S107.2pg.43",
        "id": "S107.2pg.43*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV7-2*01",
        "imgt family name": "IGHV7",
        "imgt gene name": "IGHV7-2",
        "length": 306,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTGGAGTCTGAAGGAGGCTTGGTACAGCCTGGGGGTTCTCTGAGACTCTCCTGTGCAACTTCTGGGTTCACCTTCACTGATTTCTACATGAACTGGGTCTGCCAGCCTCCAAGGAAGGCACTTGAGTGGTTGGGTTTTATTAGAAACAAAGCTAATGGTTACACAACAGACTACAGTGCATCTATGAAGGGTCGGTTCACCATCTCCAGAGATAATTCCCAAAGCATCCTCTATCTTCAAATGAACACACTGAGAACTGAGGACAGTGCCACTTATTACTGTGCAAGAGATACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2158818,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "S107.3.62*01",
        "confirmed": True,
        "end": 1917932,
        "family": "S107",
        "framed": True,
        "functionality": "F",
        "gene": "S107.3.62",
        "id": "S107.3.62*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV7-3*01",
        "imgt family name": "IGHV7",
        "imgt gene name": "IGHV7-3",
        "length": 304,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTGGAGTCTGGAGGAGGCTTGGTACAGCCTGGGGGTTCTCTGAGTCTCTCCTGTGCAGCTTCTGGATTCACCTTCACTGATTACTACATGAGCTGGGTCCGCCAGCCTCCAGGGAAGGCACTTGAGTGGTTGGGTTTTATTAGAAACAAAGCTAATGGTTACACAACAGAGTACAGTGCATCTGTGAAGGGTCGGTTCACCATCTCCAGAGATAATTCCCAAAGCATCCTCTATCTTCAAATGAATGCCCTGAGAGCTGAGGACAGTGCCACTTATTACTGTGCAAGATATA",
            "read frame": 0,
            "strand": True
        },
        "start": 1917628,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "S107.4.65*01",
        "confirmed": True,
        "end": 1848324,
        "family": "S107",
        "framed": True,
        "functionality": "F",
        "gene": "S107.4.65",
        "id": "S107.4.65*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV7-4*01",
        "imgt family name": "IGHV7",
        "imgt gene name": "IGHV7-4",
        "length": 306,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGATGGAATCTGGAGGAGGCTTGGTACAGCCTGGGGCTTCTCTGAGACTCTCCTGTGCAGCTTCTGGATTCACCTTTACTGATTACTACATGAGCTGGGTCCGCCAGCCTCCAGGGAAGGCACCTGAGTGGTTGGCTTTGATTAGAAACAAAGCTAATGGTTACACAACAGAGTATACTGCATCTGTTAAGGGTCGGTTCACCATCTCCAGAGATAATTCCCAAAACATCCTCTATCTTCAAATGAACACCCTGAGGGCTGAGGACAGTGCCACTTATTACTGTGTAAAAGCTGTA",
            "read frame": 0,
            "strand": True
        },
        "start": 1848018,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "SM7.1.44*01",
        "confirmed": True,
        "end": 2139196,
        "family": "SM7",
        "framed": True,
        "functionality": "F",
        "gene": "SM7.1.44",
        "id": "SM7.1.44*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV14-1*01",
        "imgt family name": "IGHV14",
        "imgt gene name": "IGHV14-1",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTCAGCTGCAGCAGTCTGGGGCAGAGCTTGTGAGGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTCAACATTAAAGACTACTATATGCACTGGGTGAAGCAGAGGCCTGAACAGGGCCTGGAGTGGATTGGAAGGATTGATCCTGAGGATGGTGATACTGAATATGCCCCGAAGTTCCAGGGCAAGGCCACTATGACTGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTACTACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2138902,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "SM7.2.49*01",
        "confirmed": True,
        "end": 2076643,
        "family": "SM7",
        "framed": True,
        "functionality": "F",
        "gene": "SM7.2.49",
        "id": "SM7.2.49*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV14-2*01",
        "imgt family name": "IGHV14",
        "imgt gene name": "IGHV14-2",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTCAGCTGCAGCAGTCTGGGGCAGAGCTTGTGAAGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTCAACATTAAAGACTACTATATGCACTGGGTGAAGCAGAGGACTGAACAGGGCCTGGAGTGGATTGGAAGGATTGATCCTGAGGATGGTGAAACTAAATATGCCCCGAAATTCCAGGGCAAGGCCACTATAACAGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTGCTAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2076349,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "SM7.3.54*01",
        "confirmed": True,
        "end": 2011267,
        "family": "SM7",
        "framed": True,
        "functionality": "F",
        "gene": "SM7.3.54",
        "id": "SM7.3.54*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV14-3*01",
        "imgt family name": "IGHV14",
        "imgt gene name": "IGHV14-3",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTCAGCTGCAGCAGTCTGTGGCAGAGCTTGTGAGGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTCAACATTAAAAACACCTATATGCACTGGGTGAAGCAGAGGCCTGAACAGGGCCTGGAGTGGATTGGAAGGATTGATCCTGCGAATGGTAATACTAAATATGCCCCGAAGTTCCAGGGCAAGGCCACTATAACTGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCATCTATTACTGTGCTAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2010973,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "SM7.4.63*01",
        "confirmed": True,
        "end": 1894674,
        "family": "SM7",
        "framed": True,
        "functionality": "F",
        "gene": "SM7.4.63",
        "id": "SM7.4.63*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV14-4*01",
        "imgt family name": "IGHV14",
        "imgt gene name": "IGHV14-4",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTCAGCTGCAGCAGTCTGGGGCTGAGCTTGTGAGGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTTAACATTAAAGACGACTATATGCACTGGGTGAAGCAGAGGCCTGAACAGGGCCTGGAGTGGATTGGATGGATTGATCCTGAGAATGGTGATACTGAATATGCCTCGAAGTTCCAGGGCAAGGCCACTATAACAGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTACTACA",
            "read frame": 0,
            "strand": True
        },
        "start": 1894380,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VGAM3.8-1-57*01",
        "confirmed": True,
        "end": 1977184,
        "family": "VGAM3.8",
        "framed": True,
        "functionality": "F",
        "gene": "VGAM3.8-1-57",
        "id": "VGAM3.8-1-57*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV9-1*01",
        "imgt family name": "IGHV9",
        "imgt gene name": "IGHV9-1",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATCCAGTTGGTGCAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGGGTATACCTTCACAGAATATCCAATGCACTGGGTGAAGCAGGCTCCAGGAAAGGGTTTCAAGTGGATGGGCATGATATACACCGACACTGGAGAGCCAACATATGCTGAAGAGTTCAAGGGACGGTTTGCCTTCTCTTTGGAGACCTCTGCCAGCACTGCCTATTTGCAGATCAACAACCTCAAAAATGAGGACACGGCTACATATTTCTGTGTAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1976890,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VGAM3.8-2-59*01",
        "confirmed": True,
        "end": 1962111,
        "family": "VGAM3.8",
        "framed": True,
        "functionality": "F",
        "gene": "VGAM3.8-2-59",
        "id": "VGAM3.8-2-59*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV9-2*01",
        "imgt family name": "IGHV9",
        "imgt gene name": "IGHV9-2",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATCCAGTTCGTGCAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGTGTATACCTTCACAGAATATCCAATGCACTGGGTGAAGCAGGCTCCAGGAAAGGGTTTCAAGTGGATGGGCTGGATAAACACCTACTCTGGAGAGCCAACATATGCTGACGACTTCAAGGGACGGTTTGCCTTTTCTTTGGAAACCTCTGCCAGCACTGCCTATTTGCAGATCAACAACCTCAAAAATGAGGACACGGCTACATATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1961817,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VGAM3.8-3-61*01",
        "confirmed": True,
        "end": 1930420,
        "family": "VGAM3.8",
        "framed": True,
        "functionality": "F",
        "gene": "VGAM3.8-3-61",
        "id": "VGAM3.8-3-61*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV9-3*01",
        "imgt family name": "IGHV9",
        "imgt gene name": "IGHV9-3",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATCCAGTTGGTACAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGGGTATACCTTCACAACCTATGGAATGAGCTGGGTGAAACAGGCTCCAGGAAAGGGTTTAAAGTGGATGGGCTGGATAAACACCTACTCTGGAGTGCCAACATATGCTGATGACTTCAAGGGACGGTTTGCCTTCTCTTTGGAAACCTCTGCCAGCACTGCCTATTTGCAGATCAACAACCTCAAAAATGAGGACACGGCTACATATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1930126,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VGAM3.8-4-71*01",
        "confirmed": True,
        "end": 1771148,
        "family": "VGAM3.8",
        "framed": True,
        "functionality": "F",
        "gene": "VGAM3.8-4-71",
        "id": "VGAM3.8-4-71*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV9-4*01",
        "imgt family name": "IGHV9",
        "imgt gene name": "IGHV9-4",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATCCAGTTGGTGCAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGGGTATACCTTCACAACTGCTGGAATGCAGTGGGTGCAAAAGATGCCAGGAAAGGGTTTTAAGTGGATTGGCTGGATAAACACCCACTCTGGAGAGCCAAAATATGCAGAAGACTTCAAGGGACGGTTTGCCTTCTCTTTGGAAACCTCTGCCAGCACTGCCTATTTACAGATAAGCAACCTCAAAAATGAGGACACGGCTACGTATTTCTGTGCGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1770854,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VH10.1.86*01",
        "confirmed": True,
        "end": 1592104,
        "family": "VH10",
        "framed": True,
        "functionality": "F",
        "gene": "VH10.1.86",
        "id": "VH10.1.86*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV10-1*01",
        "imgt family name": "IGHV10",
        "imgt gene name": "IGHV10-1",
        "length": 302,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGCAGCTTGTTGAGTCTGGTGGAGGATTGGTGCAGCCTAAAGGGTCATTGAAACTCTCATGTGCAGCCTCTGGATTCAGCTTCAATACCTACGCCATGAACTGGGTCCGCCAGGCTCCAGGAAAGGGTTTGGAATGGGTTGCTCGCATAAGAAGTAAAAGTAATAATTATGCAACATATTATGCCGATTCAGTGAAAGACAGATTCACCATCTCCAGAGATGATTCAGAAAGCATGCTCTATCTGCAAATGAACAACTTGAAAACTGAGGACACAGCCATGTATTACTGTGTGAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 1591802,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VH10.2pg.89*01",
        "confirmed": True,
        "end": 1575106,
        "family": "VH10",
        "framed": True,
        "functionality": "P",
        "gene": "VH10.2pg.89",
        "id": "VH10.2pg.89*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV10-2*01",
        "imgt family name": "IGHV10",
        "imgt gene name": "IGHV10-2",
        "length": 298,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGCATTGTGAGGTACAGCTTGTTGGGTCTGGTGAAGGATTGGTGCAGCCTAAAGGGTCACTTGTATTCAGCCTCTGGATTCATCTTCAATACATACACCTTGGAGTGGGTGTGTCAGACTCCAAGAAAGGGTCTGGAATGGGTTGCATGCATAAGAACAAAAAGTAATAATTATGCCACATATACTGATTCAGTGAAAGACAGATTCACCATCTCCAGAGATGATTCTCAAAGCATGGTCAACCTGCAAATGAACAATTTGAAAACTGAGGACATAGACCTGTATTAATTACAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1574808,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VH10.3.91*01",
        "confirmed": True,
        "end": 1547668,
        "family": "VH10",
        "framed": True,
        "functionality": "F",
        "gene": "VH10.3.91",
        "id": "VH10.3.91*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV10-3*01",
        "imgt family name": "IGHV10",
        "imgt gene name": "IGHV10-3",
        "length": 302,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGCAGCTTGTTGAGTCTGGTGGAGGATTGGTGCAGCCTAAAGGATCATTGAAACTCTCATGTGCCGCCTCTGGTTTCACCTTCAATACCTATGCCATGCACTGGGTCCGCCAGGCTCCAGGAAAGGGTTTGGAATGGGTTGCTCGCATAAGAAGTAAAAGTAGTAATTATGCAACATATTATGCCGATTCAGTGAAAGACAGATTCACCATCTCCAGAGATGATTCACAAAGCATGCTCTATCTGCAAATGAACAACCTGAAAACTGAGGACACAGCCATGTATTACTGTGTGAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1547366,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VH11.1.48*01",
        "confirmed": True,
        "end": 2089233,
        "family": "VH11",
        "framed": True,
        "functionality": "F",
        "gene": "VH11.1.48",
        "id": "VH11.1.48*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV11-1*01",
        "imgt family name": "IGHV11",
        "imgt gene name": "IGHV11-1",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGCAGCTGTTGGAGACTGGAGAAGGCTTGGTGCCACCTGGGGGGTCACGGGGACTCTCTTGTGAAGGCTCAGGGTTCACTTTTAGTGGCTTCTGGATGAGCTGGGTTCGACAGACACCTGGGAAGACCCTGGAGTGGATTGGAGACATTAATTCTGATGGCAGTGCAATAAACTACGCACCATCCATAAAGGATCGCTTCACTATCTTCAGAGACAATGACAAGAGCACCCTGTACCTGCAGATGAGCAATGTGCGATCTGAGGACACAGCCACGTATTTCTGTATGAGATA",
            "read frame": 0,
            "strand": True
        },
        "start": 2088937,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VH11.2.53*01",
        "confirmed": True,
        "end": 2022871,
        "family": "VH11",
        "framed": True,
        "functionality": "F",
        "gene": "VH11.2.53",
        "id": "VH11.2.53*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV11-2*01",
        "imgt family name": "IGHV11",
        "imgt gene name": "IGHV11-2",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGCAGCTGTTGGAGACTGGAGGAGGCTTGGTGCAACCTGGGGGGTCACGGGGACTCTCTTGTGAAGGCTCAGGGTTCACTTTTAGTGGCTTCTGGATGAGCTGGGTTCGACAGACACCTGGGAAGACCCTGGAGTGGATTGGAGACATTAATTCTGATGGCAGTGCAATAAACTACGCACCATCCATAAAGGATCGATTCACTATCTTCAGAGACAATGACAAGAGCACCCTGTACCTGCAGATGAGCAATGTGCGATCGGAGGACACAGCCACGTATTTCTGTATGAGATA",
            "read frame": 0,
            "strand": True
        },
        "start": 2022575,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.11.58*01",
        "confirmed": True,
        "end": 1963853,
        "family": "VH12",
        "framed": True,
        "functionality": "P",
        "gene": "PG.11.58",
        "id": "PG.11.58*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV12-1*01",
        "imgt family name": "IGHV12",
        "imgt gene name": "IGHV12-1",
        "length": 305,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "AAGATTTAGCTTAAGGAGTCAGGATCTGCTCTCATCAAGCCATCACAGCCACTGTCTCTCACCTGCACAGTCTCTGGATTCTCCATCACAAGTAGTAGTTATTGCTGGCACTGGATCTGCCAGCCCCCAGGAAAGGGGTTAGAATGGATGGGGCACATATGTTATGAAGGTTCACTAAACTATAGTCCATCCCTCAAAAGCCGCAGCACCATCTCCAGAGACACATCTCTGAACAAATTCTTTATCCAGCTGAGCTCTCTGACTGATGAGGACACAGTCATGTACTACTGTTCTAGGGAAAACCA",
            "read frame": 0,
            "strand": True
        },
        "start": 1963548,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.12.60*01",
        "confirmed": True,
        "end": 1943579,
        "family": "VH12",
        "framed": True,
        "functionality": "P",
        "gene": "PG.12.60",
        "id": "PG.12.60*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV12-2*01",
        "imgt family name": "IGHV12",
        "imgt gene name": "IGHV12-2",
        "length": 305,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATTCAGCTTAAGGAATCAGGACCTGCTGTCATCAAGCCATCATAGTCACTGTCTCTCACCTGCACAGTCTCTGGATTCTCCATCACTAGTAGTGGTTTTTGCTGGCACTGGATATGCCAGCCCCCAGGAAAGGGGTTAGAGTGGATGGGGCGCATATGTTATGAAGGTTCCATATACTATAGTCCATCCCTCAAAAGTCCCAGCACCATCTCCAGAGACATATCACTGAATAAATTCTTTATCCAGCTGAGCTCTGTAACTGATGAGGACACAGCCATGTACTACTATTCCAGGGAAAACCA",
            "read frame": 0,
            "strand": True
        },
        "start": 1943274,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VH12.1.78*01",
        "confirmed": True,
        "end": 1704589,
        "family": "VH12",
        "framed": True,
        "functionality": "F",
        "gene": "VH12.1.78",
        "id": "VH12.1.78*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV12-3*01",
        "imgt family name": "IGHV12",
        "imgt gene name": "IGHV12-3",
        "length": 300,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATGCAGCTTCAGGAGTCAGGACCTGGCCTGGTGAAACCCTCACAGTCACTCTTCCTTACCTGCTCTATTACTGGTTTCCCCATCACCAGTGGTTACTACTGGATCTGGATCCGTCAGTCACCTGGGAAACCCCTAGAATGGATGGGGTACATCACTCATAGTGGGGAAACTTTCTACAACCCATCTCTCCAGAGCCCCATCTCCATTACTAGAGAAACGTCAAAGAACCAGTTCTTCCTCCAATTGAACTCTGTGACCACAGAGGACACAGCCATGTATTACTGTGCAGGAGACAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1704289,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "IGHV15-1*01",
        "confirmed": True,
        "end": 1913087,
        "family": "VH15",
        "framed": True,
        "functionality": "P",
        "gene": "IGHV15-1",
        "id": "IGHV15-1*01*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV15-1*01",
        "imgt family name": "IGHV15",
        "imgt gene name": "IGHV15-1",
        "length": 291,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CACTCCATTGTTTAGCTGCAGCAGTCTGAAGCTGCGCTGAGGACTCTTGGAGCTTCAGAATAGGTGCCCTTCAAATGCTGTGATATGGGTAGCTTTCCCTTTTCCTTTAGGAATTGCATGACATAGAACCCTGAATGTAATGGAAACATAGACCCAAGCATGAGAAATACACTCTATGGACAGAATTAACAAGGCAGAGTCACAATGGATGCAGACAAAATGTCCAAAACTGCCTACATGTAGGTCAATAATCTGACAGCTGAGGACTCTTCTATCACTGCAGAAGGAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1912796,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VH15.1.95*01",
        "confirmed": True,
        "end": 1506530,
        "family": "VH15",
        "framed": True,
        "functionality": "F",
        "gene": "VH15.1.95",
        "id": "VH15.1.95*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV15-2*01",
        "imgt family name": "IGHV15",
        "imgt gene name": "IGHV15-2",
        "length": 297,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTCACCTACAACAGTCTGGTTCTGAACTGAGGAGTCCTGGGTCTTCAGTAAAGCTGTCATGCAAGGATTTTGATTCAGAAGTCTTCCCTATTGCTTATATGAGTTGGGTAAGGCAGAAGCCTGGGCATGGATTTGAATGGATTGGAGGCATACTCCCAAGTATTGGTAGAACAATCTATGGAGAGAAGTTTGAGGACAAAGCCACATTGGATGCAGACACACTGTCCAACACAGCCTACTTGGAGCTCAACAGTCTGACATCCGAGGACTCTGCTATCTACTACTGTGCAAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1506233,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "VH16.1.55*01",
        "confirmed": True,
        "end": 2002284,
        "family": "VH16",
        "framed": True,
        "functionality": "F",
        "gene": "VH16.1.55",
        "id": "VH16.1.55*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV16-1*01",
        "imgt family name": "IGHV16",
        "imgt gene name": "IGHV16-1",
        "length": 299,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGCAGCTGGTGGAATCTGGAGGCAGCTTGGGACAGCCTGGAGGGTCCACTAAACTCTCTTGTGAAGAAGCATCTGGATTCACTTTCAGTGATCATTGGATGGACTGGTTTCGCCAAGCCCCAGGCATGAGGCTAGAATGGTTAGCAAATACAAACCATGATGAGAGTGGAAAAGGCTATGCAGAGTCTGTGAAAGACAGATTCTCCATCTCCAGAGACAATTCTGAGAACTTATTGTATCTACAAATGAACAGTCTGAGAAACGAAGACACGGCTCTGTATTATTGTGCCAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2001985,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "X24.1pg.45*01",
        "confirmed": True,
        "end": 2122865,
        "family": "X24",
        "framed": True,
        "functionality": "F",
        "gene": "X24.1pg.45",
        "id": "X24.1pg.45*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV4-1*01",
        "imgt family name": "IGHV4",
        "imgt gene name": "IGHV4-1",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTTCTCCAGTCTGGAGGTGGCCTGGTGCAGCCTGGAGGATCCCTGAAACTCTCCTGTGCAGCCTCAGGAATCGATTTTAGTAGATACTGGATGAGTTGGGTTCGGCGGGCTCCAGGGAAAGGACTAGAATGGATTGGAGAAATTAATCCAGATAGCAGTACAATAAACTATGCACCATCTCTAAAGGATAAATTCATCATCTCCAGAGACAACGCCAAAAATACGCTGTACCTGCAAATGAGCAAAGTGAGATCTGAGGACACAGCCCTTTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2122571,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "X24.2.50*01",
        "confirmed": True,
        "end": 2057968,
        "family": "X24",
        "framed": True,
        "functionality": "P",
        "gene": "X24.2.50",
        "id": "X24.2.50*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV4-2*01",
        "imgt family name": "IGHV4",
        "imgt gene name": "IGHV4-2",
        "length": 296,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTTCTCGAGTCTGGAGGTGGCCTGGTGCAGCCTGGAGGATCCCTGAAACTCTCCTGTGCAGCCTCAGGATTCGATTTTAGTAAAGACTGGATGAGTTGGGTCCGGCAGGCTACAGGGAAAGGGCTAGAATGAATTGGAGAAATTAATCCAGGTAGCAGTACGATAAACTATACTCCATCTCTAAAGGATAAATTCATCATCTCCAGAGACAACGCCAAAAATACGCTGTACCTGCAAATGAGCAAAGTGAGATCTGAGGACACAGCCCTTTATTACTGTGCAAGACC",
            "read frame": 0,
            "strand": True
        },
        "start": 2057672,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.13.69*01",
        "end": 1786429,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.13.69",
        "id": "PG.13.69*01-C57BL/6",
        "length": 310,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGAGGTGCAGGTGGTAGAATCCGGAGGCAGCTTGATTCAGCTGGGGGGTGTGTGTCGATTAATCTCTCTTGTGAAGCTTCTGGATTCACCTTCAGTAATTACTGGATGGACTGGATTTGCCAAGCCTCAAGGAAGGGGTTAGAGTTGTTAGCAAAATCAAGAGGAATTCTATTGTCACAGCATTCTGGTGTTGCCCTGCACAAAAAGGAGAGAGATAAAGACCCCACGAAGCCTGATTGTTGAGTTTTTATACAGTTTTCAGAGCAGCACCCATTAGGAACAATGTTATTGGTAGAACGGTGTGACTTTT",
            "read frame": 0,
            "strand": True
        },
        "start": 1786119,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.19.161*01",
        "end": 658151,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.19.161",
        "id": "PG.19.161*01-C57BL/6",
        "length": 195,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGGGATGGAGCTATATCATCTTCTTCCTTGTAGCAACAGCTATATGTAAGGGTCTCACAGTAGCAGATTTGAGAGGTCTGGCCATACACTCACATGACAATGACATCCATTCTATCCTTCCATCCACAGGTATCCACTCCCAGGTATAGCTGCATCAATCTGGGGCTGAGCTGGTGAAGCCTGGGATCTCTGTTA",
            "read frame": 0,
            "strand": True
        },
        "start": 657956,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.3.23*01",
        "end": 2352251,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.3.23",
        "id": "PG.3.23*01-C57BL/6",
        "length": 314,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGGCTCAGCTGTGTTTTCCTTGTCCTCATTTTAAGAGGTAATTTGTAGAAATAAGATCCTGCCTGTTTTGTGTACAGGAGAAATAGAAATTTTTTTCTTTCCTCTACTGTGTTTTGTTTTGTTAGTGACAGTTTACAAATAAGCATTCTCTGTTGTGAGGTGTCCAGTGTGAGGTGAAGATGGTGGAGTCTGTGGGAGGCTTAGTGCAGCCTGTAGGGTCATTGAAACCCTCCTGTGCAGCCTCGGGATTCATTCTCACTGACTACTGAATGACCTGGATCCTTCAGGCTTCAAAGAAAAGGATGGAGAGGGTG",
            "read frame": 0,
            "strand": True
        },
        "start": 2351937,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.4.28*01",
        "end": 2303285,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.4.28",
        "id": "PG.4.28*01-C57BL/6",
        "length": 316,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "AAAGGTAATTCATAAAGATAAGATTCTGTCTGTTGTGTGCACATGAGAAACAGAAAAATTGTATTGTTTCTCTATTTGGTTTTGTTTTGTTAGTGACAGTTTCTGACTCAGAATTCTCTGTTTGAAGGTGCCCAGTGTGAGGTGAAGCTTGAGAAGTCTAGGGGAGGCTTAGTGCAGCCTGGAAGGTCCGTGATACGCTCATGTGCAGCCTCTTGATGCACTGCCAGTGACGACTGGTTTGCCTGGGTTTGCCAGGCTCCAAAGAAGGGGCTGCAGTGGGGGTGGGGCATAATTTTTCATGGCTGTGGTAGCCCCT",
            "read frame": 0,
            "strand": True
        },
        "start": 2302969,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.5.30*01",
        "end": 2299106,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.5.30",
        "id": "PG.5.30*01-C57BL/6",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GACCACTCACCATAGACTTTGGGCTCAGCTTTGCTTTCCTTGTTCTTATTTTAAAAGGTAATTCTTAGAAATAAGATCCTGCCTGTTTTGTGTACATGAGAAATAGAAAAATTTGTTTTCTTTCCTCTATTTTGTTTTGTTTTGTTTTGTTAGTGACAGTTTCCAAATCAGCATTCTCTGCTTTGAGGTGCCCAGTGTGAGGTGAAGCTGGTGGAGTCTGGGGGAGGCTGAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCATTCTCACTGACTACTGAATGACCT",
            "read frame": 0,
            "strand": True
        },
        "start": 2298802,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.6.32*01",
        "end": 2277088,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.6.32",
        "id": "PG.6.32*01-C57BL/6",
        "length": 297,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "AGGTGAAGCTGGAGAAGTCTAAGGGGAGGTTTAGTGCAGCCTGGAAGGTCCATGATACTCTACTGTGCAGCCTCTGGATTCACTGTCAGCGACGACTGGTTTGCCTGGGTTTGCCAGGCTCCAAAGAAGGGGCTGGAGATGGGGGGGGGGCATAATTTTTCATGGTGGTGGTAGCCCCTCTTATGCAGACACCTTGAAGAAGTGGGTTGGGCCTAACATATTCAGAATCAATATTTAAAAATTCTAATCCTTGAGTATGTCACTTTTGACCAAGAATATATGAAATATGTTACTGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 2276791,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.8.47*01",
        "end": 2102619,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.8.47",
        "id": "PG.8.47*01-C57BL/6",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TACCACTTGTAGAATCTGAAGGCAGCTTGTTATAGCCCGGAGTGTTCATTAAATTCTCTTCTGAAGCGACTGGATCCACCTTCAGTGATGCACTGGATTCGCCAAGCCCCAGGGAAGGGGCTAGAGTGGGTCAGCAGATATAAAATAATAGTGGAGTGTAAAAACTATGCAGAGTCTGTGAAGGATAGATTTGCCATCTACAGGAACAATTCTAAGAGCTTATTGTAACTACAAATAAACAAATCTAAGAAGTGAGGACACTGCCATGTGTTACTGTGTGAGAGATAAAGTGAGGAAAGTAAAT",
            "read frame": 0,
            "strand": True
        },
        "start": 2102315,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "PG.9.52*01",
        "end": 2033174,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.9.52",
        "id": "PG.9.52*01-C57BL/6",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "AGATACCACTTGTAGAATCTGAAGGCAGCTTGTTATAGCCTGGAGTGTTCATTAAATTCTCTTCTGAAGCGACTGGATCCACCTTCAGTGATGCACTGGATTCACCAAGCCCCAGGGAAGGGGCTAGAGTGGGTCAGCAGATATAAAATAATAGTGGAGTGTAAAAACTATGCAGAGTCTGTGAAGGATAGATTTGCCATCTACAGGAACAATTCTAAGAGCTTATTGTATCTACAAATAAACAAATCTAAGAAGTGAGGACACTGCCATGTGTTACTGTGTGAGAGATACAGTGAGGAAAGTAAAT",
            "read frame": 0,
            "strand": True
        },
        "start": 2032867,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "accession": "BN000872",
        "allele": "J558.20pg.110*01",
        "confirmed": True,
        "end": 1334835,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.20pg.110",
        "id": "J558.20pg.110*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-21-1*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-21-1",
        "length": 294,

        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCAGCTGCAACAGTCTGGACCTGAGTTGGTGAAGCCTGGGGCTTTAGTGAAGTTATCCTGCAAGGTTTCTGGATTCACATTCACTGACTAATACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACATGTTTATCCTTACAATGGTGGTACTAGCTACAACCAGAAATTCAAGGGCAAGGCCACATTGACTGTCGACAATACCTCCAGCACAGCCTACATGGAGCTCGGCAGCCTGACTTCTGAGGACTCTGCGGTCTATTACTCTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1334541,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    }
]

bad = [
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2106771,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2106764
            },
            "3 nonamer": {
                "end": 2106803,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2106794
            },
            "end": 2106764,
            "length": 307,
            "sequence": {
                "nucleotide": "GTATCCTGTCTGATGTGCAGCTTCAGGAGTCAGGACCTGGCATGGTGAAACCTTCTCAGTCACTTTCCCTCACCTGCACTGTCACTGGCTACTCCATCACCAGTGGTTATGACTGGCACTGGATCCGACATTTTCCAGGAAACAAACTGGAGTGGATGGGCTACATAAGCTACAGTGGTAGCACTAACTACAACCCATCCCTCAAAAGTCGAATCTCCATCACTCATGACACATCTAAGAACCATTTCTTCCTGAAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 2106457
        },
        "accession": "BN000872",
        "allele": "36-60.1.46*01",
        "confirmed": True,
        "end": 2106764,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.1.46",
        "id": "36-60.1.46*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-1*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-1",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCTGTCTGATGTGCAGCTTCAGGAGTCAGGACCTGGCATGGTGAAACCTTCTCAGTCACTTTCCCTCACCTGCACTGTCACTGGCTACTCCATCACCAGTGGTTATGACTGGCACTGGATCCGACATTTTCCAGGAAACAAACTGGAGTGGATGGGCTACATAAGCTACAGTGGTAGCACTAACTACAACCCATCCCTCAAAAGTCGAATCTCCATCACTCATGACACATCTAAGAACCATTTCTTCCTGAAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 2106457,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2037316,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2037309
            },
            "3 nonamer": {
                "end": 2037348,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2037339
            },
            "end": 2037309,
            "length": 306,
            "sequence": {
                "nucleotide": "GTATCCTGTCTGATGTGCAGCTTCAGGAGTCAGGACCTGGCCTGGTGAAACCTTCTCAGTCACTTTCCCTCACCTGCACTGTCACTGGCTACTCCATCACCAGTGGTTATGACTGGAACTGGATCCGGCAGTTTCCAGGAAACAAACTGGAGTGGATGGGCTACATAAGCTACAGTGTAGCAAGTACTACAACCCATCTCTCAAAAGTCGAATCTCTATCACTCGAGACACATCCAAGAACCAGTTCTCCCTGGAATTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 2037003
        },
        "accession": "BN000872",
        "allele": "36-60.2pg.51*01",
        "confirmed": True,
        "end": 2037309,
        "family": "36-60",
        "framed": False,
        "functionality": "P",
        "gene": "36-60.2pg.51",
        "id": "36-60.2pg.51*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-2*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-2",
        "length": 306,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCTGTCTGATGTGCAGCTTCAGGAGTCAGGACCTGGCCTGGTGAAACCTTCTCAGTCACTTTCCCTCACCTGCACTGTCACTGGCTACTCCATCACCAGTGGTTATGACTGGAACTGGATCCGGCAGTTTCCAGGAAACAAACTGGAGTGGATGGGCTACATAAGCTACAGTGTAGCAAGTACTACAACCCATCTCTCAAAAGTCGAATCTCTATCACTCGAGACACATCCAAGAACCAGTTCTCCCTGGAATTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2037003,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1874681,
                "length": 8,
                "sequence": {
                    "nucleotide": "CACAATGT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1874673
            },
            "3 nonamer": {
                "end": 1874712,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACATAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1874703
            },
            "end": 1874673,
            "length": 307,
            "sequence": {
                "nucleotide": "GTATCCTGTCTGATGTGCAGCTTCAGGAGTCAGGACCTAGCCTGGTGAGACCTTCTCAGACACTCTCCCTTACCTGCACTGTTACTGGCTTCTCCATCAACAGTGATTGTTACTGGATCTGGATCCGGCAGTTTCCAGGAAACAAACTGGAGTACATCGGGTACACATTCTACAGTGGTATCACTTACTACAACCCATCTCTTGAAAGTCGAACGTACATAACGCGTGACACATCTAAGAACCAGTTCTCACTGAAGTTGAGTTCTGTGACTACTGAGGACACAGCCACTTACTACTGTGCGAGAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1874366
        },
        "accession": "BN000872",
        "allele": "36-60.3.64*01",
        "confirmed": True,
        "end": 1874673,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.3.64",
        "id": "36-60.3.64*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-3*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-3",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCTGTCTGATGTGCAGCTTCAGGAGTCAGGACCTAGCCTGGTGAGACCTTCTCAGACACTCTCCCTTACCTGCACTGTTACTGGCTTCTCCATCAACAGTGATTGTTACTGGATCTGGATCCGGCAGTTTCCAGGAAACAAACTGGAGTACATCGGGTACACATTCTACAGTGGTATCACTTACTACAACCCATCTCTTGAAAGTCGAACGTACATAACGCGTGACACATCTAAGAACCAGTTCTCACTGAAGTTGAGTTCTGTGACTACTGAGGACACAGCCACTTACTACTGTGCGAGAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1874366,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1817502,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTA",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1817495
            },
            "3 nonamer": {
                "end": 1817534,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1817525
            },
            "end": 1817495,
            "length": 310,
            "sequence": {
                "nucleotide": "GTATCTTGTCTGATGTACAGCTTCAGGAGTCAGGACCTGCCCTGGTGAAGCCTTCTCAGACAGTGTCCCTCACCTGCACTGTCACTGGCTACTCTATCACTAATGGTAATCACTGGTGGAACTGGATCCGGCAGGTTTCAGGAAGCAAACTGGAGTGGATAGGGTACATAAGCTCCAGTGGTAGCACTGACAGCAATCCATCTCTCAAAAGTCGAATCTCCATCACTAGAGACACTTCCAAGAACCAGTTATTCCTGCAGTTGAACTCTGTGACTACTGAAGATATAGCCACATATTACTGTGCAAGAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1817185
        },
        "accession": "BN000872",
        "allele": "36-60.4.66*01",
        "confirmed": True,
        "end": 1817495,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.4.66",
        "id": "36-60.4.66*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-4*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-4",
        "length": 310,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCTTGTCTGATGTACAGCTTCAGGAGTCAGGACCTGCCCTGGTGAAGCCTTCTCAGACAGTGTCCCTCACCTGCACTGTCACTGGCTACTCTATCACTAATGGTAATCACTGGTGGAACTGGATCCGGCAGGTTTCAGGAAGCAAACTGGAGTGGATAGGGTACATAAGCTCCAGTGGTAGCACTGACAGCAATCCATCTCTCAAAAGTCGAATCTCCATCACTAGAGACACTTCCAAGAACCAGTTATTCCTGCAGTTGAACTCTGTGACTACTGAAGATATAGCCACATATTACTGTGCAAGAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1817185,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1808467,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1808460
            },
            "3 nonamer": {
                "end": 1808499,
                "length": 9,
                "sequence": {
                    "nucleotide": "ATGTAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1808490
            },
            "end": 1808460,
            "length": 310,
            "sequence": {
                "nucleotide": "GTATCCTGTCTGATGTGCAGCTTCAGGAGTCAGGACCTGGCCTGGTGAAACCTTCTCAGACAGTGTTCCTCACCTGCACTGTCACTGGCATTTCCATCACCACTGGAAATTACAGGTGGAGCTGGATCCGGCAGTTTCCAGGAAACAAACTGGAGTGGATAGGGTACATATACTACAGTGGTACCATTACCTACAATCCATCTCTCACAAGTCGAACCACCATCACTAGAGACACTCCCAAGAACCAGTTCTTCCTGGAAATGAACTCTTTGACTGCTGAGGACACAGCCACATACTACTGTGCACGAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1808150
        },
        "accession": "BN000872",
        "allele": "36-60.5.67*01",
        "confirmed": True,
        "end": 1808460,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.5.67",
        "id": "36-60.5.67*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-5*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-5",
        "length": 310,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCTGTCTGATGTGCAGCTTCAGGAGTCAGGACCTGGCCTGGTGAAACCTTCTCAGACAGTGTTCCTCACCTGCACTGTCACTGGCATTTCCATCACCACTGGAAATTACAGGTGGAGCTGGATCCGGCAGTTTCCAGGAAACAAACTGGAGTGGATAGGGTACATATACTACAGTGGTACCATTACCTACAATCCATCTCTCACAAGTCGAACCACCATCACTAGAGACACTCCCAAGAACCAGTTCTTCCTGGAAATGAACTCTTTGACTGCTGAGGACACAGCCACATACTACTGTGCACGAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1808150,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1782965,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1782958
            },
            "3 nonamer": {
                "end": 1782997,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACATAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1782988
            },
            "end": 1782958,
            "length": 307,
            "sequence": {
                "nucleotide": "GTATCCTGTCTGATGTACAGCTTCAGGAGTCAGGACCTGGCCTCGTGAAACCTTCTCAGTCTCTGTCTCTCACCTGCTCTGTCACTGGCTACTCCATCACCAGTGGTTATTACTGGAACTGGATCCGGCAGTTTCCAGGAAACAAACTGGAATGGATGGGCTACATAAGCTACGATGGTAGCAATAACTACAACCCATCTCTCAAAAATCGAATCTCCATCACTCGTGACACATCTAAGAACCAGTTTTTCCTGAAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1782651
        },
        "accession": "BN000872",
        "allele": "36-60.6.70*01",
        "confirmed": True,
        "end": 1782958,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.6.70",
        "id": "36-60.6.70*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-6*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-6",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCTGTCTGATGTACAGCTTCAGGAGTCAGGACCTGGCCTCGTGAAACCTTCTCAGTCTCTGTCTCTCACCTGCTCTGTCACTGGCTACTCCATCACCAGTGGTTATTACTGGAACTGGATCCGGCAGTTTCCAGGAAACAAACTGGAATGGATGGGCTACATAAGCTACGATGGTAGCAATAACTACAACCCATCTCTCAAAAATCGAATCTCCATCACTCGTGACACATCTAAGAACCAGTTTTTCCTGAAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1782651,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1756932,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAAGG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1756925
            },
            "3 nonamer": {
                "end": 1756964,
                "length": 9,
                "sequence": {
                    "nucleotide": "AAATAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1756955
            },
            "end": 1756925,
            "length": 307,
            "sequence": {
                "nucleotide": "GTATCCTATCTGAGGTACAGCTCCAGGAGTCGAGACCTAGCCTGGTGAAACCTTCTCAGTCATTGTCCTTCACCTGCTCTGTCACTGGCTACTCCATCACCAGTGGTTATGACTGGAACTGGATCCGGCAGTTTCCAGGAAACAAACTGTAATGGATGAGGTACATACACAAGAGTGGTAGCACTAACTACAATCCATCTCTCAAAAGTCAAGTCTCCATCACTCGAGACACTTCCAAGAACCAGTTCTACCTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCATATATTACTGTGCAAACGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1756618
        },
        "accession": "BN000872",
        "allele": "36-60.7pg.72*01",
        "confirmed": True,
        "end": 1756925,
        "family": "36-60",
        "framed": True,
        "functionality": "P",
        "gene": "36-60.7pg.72",
        "id": "36-60.7pg.72*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-7*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-7",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCTATCTGAGGTACAGCTCCAGGAGTCGAGACCTAGCCTGGTGAAACCTTCTCAGTCATTGTCCTTCACCTGCTCTGTCACTGGCTACTCCATCACCAGTGGTTATGACTGGAACTGGATCCGGCAGTTTCCAGGAAACAAACTGTAATGGATGAGGTACATACACAAGAGTGGTAGCACTAACTACAATCCATCTCTCAAAAGTCAAGTCTCCATCACTCGAGACACTTCCAAGAACCAGTTCTACCTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCATATATTACTGTGCAAACGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1756618,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1748742,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1748735
            },
            "3 nonamer": {
                "end": 1748774,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACATAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1748765
            },
            "end": 1748735,
            "length": 304,
            "sequence": {
                "nucleotide": "GTATCCTGTCAGAGGTGCAGCTTCAGGAGTCAGGACCTGGCCTGGCAAAACCTTCTCAGACTCTGTCCCTCACCTGTTCTGTCACTGGCTACTCCATCACCAGTGATTACTGGAACTGGATCCGGAAATTCCCAGGGAATAAACTTGAGTACATGGGGTACATAAGCTACAGTGGTAGCACTTACTACAATCCATCTCTCAAAAGTCGAATCTCCATAACTCGAGACACATCCAAGAACCAGTATTACCTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGATA",
                "read frame": 2,
                "strand": True
            },
            "start": 1748431
        },
        "accession": "BN000872",
        "allele": "36-60.8.74*01",
        "confirmed": True,
        "end": 1748735,
        "family": "36-60",
        "framed": True,
        "functionality": "F",
        "gene": "36-60.8.74",
        "id": "36-60.8.74*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV3-8*01",
        "imgt family name": "IGHV3",
        "imgt gene name": "IGHV3-8",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCTGTCAGAGGTGCAGCTTCAGGAGTCAGGACCTGGCCTGGCAAAACCTTCTCAGACTCTGTCCCTCACCTGTTCTGTCACTGGCTACTCCATCACCAGTGATTACTGGAACTGGATCCGGAAATTCCCAGGGAATAAACTTGAGTACATGGGGTACATAAGCTACAGTGGTAGCACTTACTACAATCCATCTCTCAAAAGTCGAATCTCCATAACTCGAGACACATCCAAGAACCAGTATTACCTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGATA",
            "read frame": 2,
            "strand": True
        },
        "start": 1748431,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1803516,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1803509
            },
            "3 nonamer": {
                "end": 1803552,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1803543
            },
            "end": 1803509,
            "length": 311,
            "sequence": {
                "nucleotide": "GTATTCACTGTGAGGTACAGCTGGTAGAGACAGGAGGAGGCTTGGTGCAGCCTGGAAACTCTCTAAAACTTTCCTGTGCCACTTCGGGATACCCTTTTTATGACTACTGGATGGATTGGGTCCGCCATTCTCCAGAAAAGGGGCTGGAGTGGGTTGCTCGAATTGCAACAAAAACTCATAATTATGCAACGTACTATGCAGAGTCTGTGAAAGGCCGATTCATCGTCTCAAGAGATGATTCCAAAAGCAGTGCATACATGCAGATGAACAGCTTAAGAAAGGAAGACACTGCCGTTTATTACTGTGCAAGAGAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1803194
        },
        "accession": "BN000872",
        "allele": "3609N.1pg.68*01",
        "confirmed": True,
        "end": 1803509,
        "family": "3609N",
        "framed": True,
        "functionality": "P",
        "gene": "3609N.1pg.68",
        "id": "3609N.1pg.68*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV13-1*01",
        "imgt family name": "IGHV13",
        "imgt gene name": "IGHV13-1",
        "length": 315,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATTCACTGTGAGGTACAGCTGGTAGAGACAGGAGGAGGCTTGGTGCAGCCTGGAAACTCTCTAAAACTTTCCTGTGCCACTTCGGGATACCCTTTTTATGACTACTGGATGGATTGGGTCCGCCATTCTCCAGAAAAGGGGCTGGAGTGGGTTGCTCGAATTGCAACAAAAACTCATAATTATGCAACGTACTATGCAGAGTCTGTGAAAGGCCGATTCATCGTCTCAAGAGATGATTCCAAAAGCAGTGCATACATGCAGATGAACAGCTTAAGAAAGGAAGACACTGCCGTTTATTACTGTGCAAGAGAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1803194,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1713355,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACACTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1713348
            },
            "3 nonamer": {
                "end": 1713387,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1713378
            },
            "end": 1713348,
            "length": 311,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTCAGGTGCAGCTTGTAGAGACCGGGGGAGGCTTGGTGAGGCCTGGAAATTCTCTGAAACTCTCCTGTGTTACCTCGGGATTCACTTTCAGTAACTACCGGATGCACTGGCTTCGCCAGCCTCCAGGGAAGAGGCTGGAGTGGATTGCTGTAATTACAGTCAAATCTGATAATTATGGAGCAAATTATGCAGAGTCTGTGAAAGGCAGATTCGCCATTTCAAGAGATGATTCAAAAAGCAGTGTCTACCTAGAGATGAACAGATTAAGAGAGGAAGACACTGCCACTTATTTTTGTAGTAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1713037
        },
        "accession": "BN000872",
        "allele": "3609N.2.77*01",
        "confirmed": True,
        "end": 1713348,
        "family": "3609N",
        "framed": True,
        "functionality": "F",
        "gene": "3609N.2.77",
        "id": "3609N.2.77*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV13-2*01",
        "imgt family name": "IGHV13",
        "imgt gene name": "IGHV13-2",
        "length": 311,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTCAGGTGCAGCTTGTAGAGACCGGGGGAGGCTTGGTGAGGCCTGGAAATTCTCTGAAACTCTCCTGTGTTACCTCGGGATTCACTTTCAGTAACTACCGGATGCACTGGCTTCGCCAGCCTCCAGGGAAGAGGCTGGAGTGGATTGCTGTAATTACAGTCAAATCTGATAATTATGGAGCAAATTATGCAGAGTCTGTGAAAGGCAGATTCGCCATTTCAAGAGATGATTCAAAAAGCAGTGTCTACCTAGAGATGAACAGATTAAGAGAGGAAGACACTGCCACTTATTTTTGTAGTAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1713037,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1608826,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1608819
            },
            "3 nonamer": {
                "end": 1608857,
                "length": 9,
                "sequence": {
                    "nucleotide": "TGCAATATT",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1608848
            },
            "end": 1608819,
            "length": 297,
            "sequence": {
                "nucleotide": "GTTACTCTGAAAGTGTCTGGCCCTGGGATATTGCAGCCATCACAGACTCTCAGCCTGGCCTGTACTTTCTCTGGGATTTCACTGAGTACTTCTGGTATGGGTTTGAGCTGGCTTCGTAAGCCCTCAGGGAAGGCTTTAGAGTGGCTGGCAAGCATTTGGAATAATGATAACTACTACAACCCATCTTTGAAGAGCCGGCTCACAATCTCCAAGGAGACCTCCAACTACCAAGTATTCCTTAAACTCACCAGTGTGGACACTGCAGATTCTGCCACATACTACGGTGCTTGGAGAGAG",
                "read frame": 0,
                "strand": True
            },
            "start": 1608522
        },
        "accession": "BN000872",
        "allele": "3609.1.84*01",
        "confirmed": True,
        "end": 1608819,
        "family": "3609v",
        "framed": True,
        "functionality": "O",
        "gene": "3609.1.84",
        "id": "3609.1.84*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-2*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-2",
        "length": 297,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTTACTCTGAAAGTGTCTGGCCCTGGGATATTGCAGCCATCACAGACTCTCAGCCTGGCCTGTACTTTCTCTGGGATTTCACTGAGTACTTCTGGTATGGGTTTGAGCTGGCTTCGTAAGCCCTCAGGGAAGGCTTTAGAGTGGCTGGCAAGCATTTGGAATAATGATAACTACTACAACCCATCTTTGAAGAGCCGGCTCACAATCTCCAAGGAGACCTCCAACTACCAAGTATTCCTTAAACTCACCAGTGTGGACACTGCAGATTCTGCCACATACTACGGTGCTTGGAGAGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 1608522,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 551738,
            "length": 312,
            "sequence": {
                "nucleotide": "ATGTCCTGTCCCAAGTTACGCTTTAACAGTCTGGCCCTGGGAAATTGCAGCCCTCCCCGACCTTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTATGGAATGATGGTGAGCTGGATGTGTCAGCCTTCAGGGAAGGGTGTGGAGTGGCTGGCACACATTTGGTGGAATGATGATAAGGGCTACAACCCATCCCTGAGGAGCCGGCTGACAATCTCCAAAGATACCCCCAACACCCAGGTATTCCTTAAGATCTCCAGTGTGGTCACTGCAGATACTACCACATACTACTGTGCTTGAGGAG",
                "read frame": 1,
                "strand": True
            },
            "start": 551426
        },
        "accession": "BN000872",
        "allele": "3609.10pg.167*01",
        "confirmed": True,
        "end": 551738,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.10pg.167",
        "id": "3609.10pg.167*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-10*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-10",
        "length": 301,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAAGTTACGCTTTAACAGTCTGGCCCTGGGAAATTGCAGCCCTCCCCGACCTTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTATGGAATGATGGTGAGCTGGATGTGTCAGCCTTCAGGGAAGGGTGTGGAGTGGCTGGCACACATTTGGTGGAATGATGATAAGGGCTACAACCCATCCCTGAGGAGCCGGCTGACAATCTCCAAAGATACCCCCAACACCCAGGTATTCCTTAAGATCTCCAGTGTGGTCACTGCAGATACTACCACATACTACTGTGCTTGAGGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 551437,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 503936,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACATTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 503929
            },
            "3 nonamer": {
                "end": 503968,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAGTAATT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 503959
            },
            "end": 503929,
            "length": 312,
            "sequence": {
                "nucleotide": "ATGTCCTCTCCCAGATTACTCAGAAAGAGTCTGGCCCTGGGATATTGCAGCCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTCTGGTATGGGTGTAGGCTGGATTCATCAGCCTTCAGGGAATGGTCTGGAGTGGCTGGCACACATTTGGTGGAATGATAATAAGTACTATAACACAGCCCTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCGCCAGTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAATAG",
                "read frame": 2,
                "strand": True
            },
            "start": 503617
        },
        "accession": "BN000872",
        "allele": "3609.11.169*01",
        "confirmed": True,
        "end": 503929,
        "family": "3609v",
        "framed": True,
        "functionality": "F",
        "gene": "3609.11.169",
        "id": "3609.11.169*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-11*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-11",
        "length": 312,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATGTCCTCTCCCAGATTACTCAGAAAGAGTCTGGCCCTGGGATATTGCAGCCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTCTGGTATGGGTGTAGGCTGGATTCATCAGCCTTCAGGGAATGGTCTGGAGTGGCTGGCACACATTTGGTGGAATGATAATAAGTACTATAACACAGCCCTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCGCCAGTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAATAG",
            "read frame": 2,
            "strand": True
        },
        "start": 503617,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 423140,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACATTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 423133
            },
            "3 nonamer": {
                "end": 423172,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAGTATTT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 423163
            },
            "end": 423133,
            "length": 312,
            "sequence": {
                "nucleotide": "ATGTCCTGTCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGGATATTGCAGTCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTCTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGAAAGGGTCTGGAGTGGCTGGCACACATTTACTGGGATGATGACAAGCGCTATAACCCATCCCTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAGAAACCAGGTATTCCTCAAGATCACCAGTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAAGAG",
                "read frame": 2,
                "strand": True
            },
            "start": 422821
        },
        "accession": "BN000872",
        "allele": "3609.12.174*01",
        "confirmed": True,
        "end": 423133,
        "family": "3609v",
        "framed": True,
        "functionality": "F",
        "gene": "3609.12.174",
        "id": "3609.12.174*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-12*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-12",
        "length": 312,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATGTCCTGTCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGGATATTGCAGTCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTCTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGAAAGGGTCTGGAGTGGCTGGCACACATTTACTGGGATGATGACAAGCGCTATAACCCATCCCTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAGAAACCAGGTATTCCTCAAGATCACCAGTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAAGAG",
            "read frame": 2,
            "strand": True
        },
        "start": 422821,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 305743,
            "length": 312,
            "sequence": {
                "nucleotide": "ATGTCCTATCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGTATATTGCAGCCCTCCCAGACCCTCAGTCTGACCTGTTCTTTCTCTGGGTTTTCACTGAGCACTTTTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTGGAGTGGCTGGCACACATTTATTGGGATGATGACAAGCACTATAACCCATCCTTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCACTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAAGAG",
                "read frame": 2,
                "strand": True
            },
            "start": 305431
        },
        "accession": "BN000872",
        "allele": "3609.13pg.178*01",
        "confirmed": True,
        "end": 305743,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.13pg.178",
        "id": "3609.13pg.178*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-13*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-13",
        "length": 301,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTTACTCTGAAAGAGTCTGGCCCTGGTATATTGCAGCCCTCCCAGACCCTCAGTCTGACCTGTTCTTTCTCTGGGTTTTCACTGAGCACTTTTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTGGAGTGGCTGGCACACATTTATTGGGATGATGACAAGCACTATAACCCATCCTTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCACTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAAGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 305442,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 263628,
            "length": 312,
            "sequence": {
                "nucleotide": "ATGTCCTGTCCCAAGTTACTCTAAAAGAGTATGGCCCTGGGAAACTGTAGCCCTCTCAGACCTTCAGTCTGACTTGTACTTTCTCTGGGTTTTCACTGAGCACTTATGGAATGATGGTGAGCTGCATGTGTCAGCCTTCAGGGAAGGGTCTGGTGTGGCTGGCACTCATTTGGTGGAATAATGATAAGGGCTACAACCCATTCCTGAGGAGCCAGCTGACAATCTCCAAGGATACATCCAACAACCAGGTATTCCTCAAGATCACCAGTGTGGACCCTGCAGATACTGCCACATACTACTGTGCTTGAGGAG",
                "read frame": 0,
                "strand": True
            },
            "start": 263316
        },
        "accession": "BN000872",
        "allele": "3609.14pg.181*01",
        "confirmed": True,
        "end": 263628,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.14pg.181",
        "id": "3609.14pg.181*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-14*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-14",
        "length": 301,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAAGTTACTCTAAAAGAGTATGGCCCTGGGAAACTGTAGCCCTCTCAGACCTTCAGTCTGACTTGTACTTTCTCTGGGTTTTCACTGAGCACTTATGGAATGATGGTGAGCTGCATGTGTCAGCCTTCAGGGAAGGGTCTGGTGTGGCTGGCACTCATTTGGTGGAATAATGATAAGGGCTACAACCCATTCCTGAGGAGCCAGCTGACAATCTCCAAGGATACATCCAACAACCAGGTATTCCTCAAGATCACCAGTGTGGACCCTGCAGATACTGCCACATACTACTGTGCTTGAGGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 263327,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 228044,
            "length": 310,
            "sequence": {
                "nucleotide": "ATGTCTTGTCAGAGTTTACTCTAAAAGAGTCTGGCCCTGGGATACTGCAGACCTCCAAGACGCTCAGTCTTACTTGCTCTTCCTTGGGTTTTCACAGAGAACTTCTGGTTTGGTCATGAGCGGGATCTGTCAGCCATCAGGGAACTGTCTGAGGTAGTTGGCATTAGTTGATGGAATGATGACAAGTACTACAACCCTTCTCTGAAGAGCTGGATCACAATATCCAAGGACACCTCCACCAACCAACCATTCCTCAAGACCACCACTGTGGACACTGCAGACACTGCCACATATTAATGTGCTCCAGGAG",
                "read frame": 2,
                "strand": True
            },
            "start": 227734
        },
        "accession": "BN000872",
        "allele": "3609.15pg.183*01",
        "confirmed": True,
        "end": 228044,
        "family": "3609v",
        "framed": False,
        "functionality": "P",
        "gene": "3609.15pg.183",
        "id": "3609.15pg.183*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-15*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-15",
        "length": 299,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGTTTACTCTAAAAGAGTCTGGCCCTGGGATACTGCAGACCTCCAAGACGCTCAGTCTTACTTGCTCTTCCTTGGGTTTTCACAGAGAACTTCTGGTTTGGTCATGAGCGGGATCTGTCAGCCATCAGGGAACTGTCTGAGGTAGTTGGCATTAGTTGATGGAATGATGACAAGTACTACAACCCTTCTCTGAAGAGCTGGATCACAATATCCAAGGACACCTCCACCAACCAACCATTCCTCAAGACCACCACTGTGGACACTGCAGACACTGCCACATATTAATGTGCTCCAGGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 227745,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 214051,
            "length": 310,
            "sequence": {
                "nucleotide": "ATGTCTTGTCAGAGTTTACTCTAAAAGAGTCTGGCCCTGGGATACTGCAGACCTCCAAGACGCTCAGTCTTACTTGCTCTTCCTTGGGTTTTCACAGAGAACTTCTGGTTTGGTCATGAGTGGGATCTGTCAGCCATCAGGGAACTGTCTGAGGTAGTTGGCATTAATTGATGGAATGATGATAAGTACTACAACCCTTCTCTGAAGAGCTGGATCACAATATCCAAGGACACCTCCACCAACCAACCATTCCTCAAGACCACCACTGTGGACACTGCAGATAGTGCCACATATTAATGTGCTCCAGGAG",
                "read frame": 0,
                "strand": True
            },
            "start": 213741
        },
        "accession": "BN000872",
        "allele": "3609.16pg.185*01",
        "confirmed": True,
        "end": 214051,
        "family": "3609v",
        "framed": False,
        "functionality": "P",
        "gene": "3609.16pg.185",
        "id": "3609.16pg.185*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-16*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-16",
        "length": 299,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGTTTACTCTAAAAGAGTCTGGCCCTGGGATACTGCAGACCTCCAAGACGCTCAGTCTTACTTGCTCTTCCTTGGGTTTTCACAGAGAACTTCTGGTTTGGTCATGAGTGGGATCTGTCAGCCATCAGGGAACTGTCTGAGGTAGTTGGCATTAATTGATGGAATGATGATAAGTACTACAACCCTTCTCTGAAGAGCTGGATCACAATATCCAAGGACACCTCCACCAACCAACCATTCCTCAAGACCACCACTGTGGACACTGCAGATAGTGCCACATATTAATGTGCTCCAGGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 213752,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 1067651,
            "length": 272,
            "sequence": {
                "nucleotide": "GTTACTTTGAAAGTGTCTGGCCCTGGCGTATTGCATCCCTCCCAGGCTCACAGTCTGACTTGCTCATTGTGTTTTCACTGAGCACTTCTGGTAGGACCATGATTCATCAGCCTTCAGGGAAAGGTTTGACAACAATTTGTTGGAATTATGATAAGTGCTACAACACATCTATTAAGAGTCAGCTCACAGTCTCCAAGGACACCTCAAACAACCAAGCATTACTCCAGATCATCAATGTGGACACTGCAAATACTGACACATACTACGGTGCC",
                "read frame": 0,
                "strand": True
            },
            "start": 1067379
        },
        "accession": "BN000872",
        "allele": "3609.2pg.138*01",
        "confirmed": True,
        "end": 1067659,
        "family": "3609v",
        "framed": False,
        "functionality": "P",
        "gene": "3609.2pg.138",
        "id": "3609.2pg.138*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-3*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-3",
        "length": 295,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTCTGTCTTAGGTTACTTTGAAAGTGTCTGGCCCTGGCGTATTGCATCCCTCCCAGGCTCACAGTCTGACTTGCTCATTGTGTTTTCACTGAGCACTTCTGGTAGGACCATGATTCATCAGCCTTCAGGGAAAGGTTTGACAACAATTTGTTGGAATTATGATAAGTGCTACAACACATCTATTAAGAGTCAGCTCACAGTCTCCAAGGACACCTCAAACAACCAAGCATTACTCCAGATCATCAATGTGGACACTGCAAATACTGACACATACTACGGTGCCTCAGTTTT",
            "read frame": 0,
            "strand": True
        },
        "start": 1067364,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1047058,
                "length": 7,
                "sequence": {
                    "nucleotide": "CAGAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1047051
            },
            "3 nonamer": {
                "end": 1047090,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGCTGTGC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1047081
            },
            "end": 1047051,
            "length": 305,
            "sequence": {
                "nucleotide": "ATGTTCTGTCCCAGATTACTCTGAAACAGTCTGGCCCTGGGATAGTGCAGCCATCCCAGCCCGTCAGACTTACTTGCACTTTCTCTGGGTTTTCACTGAGCACTTCTGGTATAGGTGTAACCTGGATTCGTCAGCCCTCAGGGAAGGGTCTGGAGTGGCTGGCAACCATTTGGTGGGATGATGATAACCGCTACAACCCATCTCTAAAGAGCAGGCTCGCAGTCTCCAAAGACACCTCCAACAACCAAGCATTCCTGAATATCATCACTGTGGAAACTGCAGATACTGCCATATACTACTGTGCT",
                "read frame": 2,
                "strand": True
            },
            "start": 1046746
        },
        "accession": "BN000872",
        "allele": "3609.3.139*01",
        "confirmed": True,
        "end": 1047051,
        "family": "3609v",
        "framed": True,
        "functionality": "O",
        "gene": "3609.3.139",
        "id": "3609.3.139*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-4*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-4",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATGTTCTGTCCCAGATTACTCTGAAACAGTCTGGCCCTGGGATAGTGCAGCCATCCCAGCCCGTCAGACTTACTTGCACTTTCTCTGGGTTTTCACTGAGCACTTCTGGTATAGGTGTAACCTGGATTCGTCAGCCCTCAGGGAAGGGTCTGGAGTGGCTGGCAACCATTTGGTGGGATGATGATAACCGCTACAACCCATCTCTAAAGAGCAGGCTCGCAGTCTCCAAAGACACCTCCAACAACCAAGCATTCCTGAATATCATCACTGTGGAAACTGCAGATACTGCCATATACTACTGTGCT",
            "read frame": 2,
            "strand": True
        },
        "start": 1046746,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1003525,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACATTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1003518
            },
            "3 nonamer": {
                "end": 1003557,
                "length": 9,
                "sequence": {
                    "nucleotide": "ATAGTAATT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1003548
            },
            "end": 1003518,
            "length": 312,
            "sequence": {
                "nucleotide": "ATGTCCTGTCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGGATATTGCAGCCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTCTAATATGGGTATAGGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTAGAGTGGCTGGCACACATTTGGTGGAATGATGATAAGTACTATAACCCATCCCTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCAGTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCAAATAG",
                "read frame": 0,
                "strand": True
            },
            "start": 1003206
        },
        "accession": "BN000872",
        "allele": "3609.4.142*01",
        "confirmed": True,
        "end": 1003518,
        "family": "3609v",
        "framed": True,
        "functionality": "F",
        "gene": "3609.4.142",
        "id": "3609.4.142*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-5*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-5",
        "length": 312,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATGTCCTGTCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGGATATTGCAGCCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTCTAATATGGGTATAGGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTAGAGTGGCTGGCACACATTTGGTGGAATGATGATAAGTACTATAACCCATCCCTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCAGTGTGGACACTGCAGATACTGCCACATACTACTGTGCTCAAATAG",
            "read frame": 0,
            "strand": True
        },
        "start": 1003206,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 905308,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACGTTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 905301
            },
            "3 nonamer": {
                "end": 905340,
                "length": 9,
                "sequence": {
                    "nucleotide": "AAAGTATTT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 905331
            },
            "end": 905301,
            "length": 312,
            "sequence": {
                "nucleotide": "ATGTCCTATCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGTATATTGCAGCCCTCGCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGTACTTTTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGGAAGGATCTGGAGTGGCTGGCACACATTTATTGGGATGATGACAAGCACTATAACCCATCCTTGAAGAGCCAGCTCAGAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCACTGTGGACACTGTAGATACTGCCACATACTACTGTGCTCGAAGAG",
                "read frame": 0,
                "strand": True
            },
            "start": 904989
        },
        "accession": "BN000872",
        "allele": "3609.5.147*01",
        "confirmed": True,
        "end": 905301,
        "family": "3609v",
        "framed": True,
        "functionality": "F",
        "gene": "3609.5.147",
        "id": "3609.5.147*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-6*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-6",
        "length": 312,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATGTCCTATCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGTATATTGCAGCCCTCGCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGTACTTTTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGGAAGGATCTGGAGTGGCTGGCACACATTTATTGGGATGATGACAAGCACTATAACCCATCCTTGAAGAGCCAGCTCAGAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCACTGTGGACACTGTAGATACTGCCACATACTACTGTGCTCGAAGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 904989,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 812953,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACATTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 812946
            },
            "3 nonamer": {
                "end": 812985,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACATTATTT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 812976
            },
            "end": 812946,
            "length": 293,
            "sequence": {
                "nucleotide": "GTTATTCTGAAAGAGTCTGGCCCTGGAATATTGCAGCCGTCCCAGACCCTAAGTCTGACTTGTTCTTTCTCTGGGTTTTCATTTAGGACTTAAGGTATGGCAGTGAGCTGGATGCGTCAGCCTTTAGGGAAGGGTTTGGAGTGGCTGGCACAAATTGGGTCAGATGATAGTAAGCTCTATAACCCATCTCTGAAGAGTCGAACCACAATCTCCAAGGATACCTCCAACAACCATGTTTTCCTCAAGATCACCAGTGAGGACACTGAAGATTCTGCCACATACTACTGTGCTCACAGAG",
                "read frame": 1,
                "strand": True
            },
            "start": 812648
        },
        "accession": "BN000872",
        "allele": "3609.6pg.151*01",
        "confirmed": True,
        "end": 812946,
        "family": "3609v",
        "framed": True,
        "functionality": "P",
        "gene": "3609.6pg.151",
        "id": "3609.6pg.151*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-7*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-7",
        "length": 298,
        "organism name": "mus musculus",
        "read frame": 2,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTTATTCTGAAAGAGTCTGGCCCTGGAATATTGCAGCCGTCCCAGACCCTAAGTCTGACTTGTTCTTTCTCTGGGTTTTCATTTAGGACTTAAGGTATGGCAGTGAGCTGGATGCGTCAGCCTTTAGGGAAGGGTTTGGAGTGGCTGGCACAAATTGGGTCAGATGATAGTAAGCTCTATAACCCATCTCTGAAGAGTCGAACCACAATCTCCAAGGATACCTCCAACAACCATGTTTTCCTCAAGATCACCAGTGAGGACACTGAAGATTCTGCCACATACTACTGTGCTCACAGAG",
            "read frame": 1,
            "strand": True
        },
        "start": 812648,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 777023,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACATTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 777016
            },
            "3 nonamer": {
                "end": 777055,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAGTAATT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 777046
            },
            "end": 777016,
            "length": 312,
            "sequence": {
                "nucleotide": "ATGTCCTGTCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGGATATTGCAGCCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTTTGGTATGGGTGTAGGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTGGAGTGGCTGGCACACATTTGGTGGGATGATGATAAGTACTATAACCCAGCCCTGAAGAGTCGGCTCACAATCTCCAAGGATACCTCCAAAAACCAGGTATTCCTCAAGATCGCCAATGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAATAG",
                "read frame": 2,
                "strand": True
            },
            "start": 776704
        },
        "accession": "BN000872",
        "allele": "3609.7.153*01",
        "confirmed": True,
        "end": 777016,
        "family": "3609v",
        "framed": True,
        "functionality": "F",
        "gene": "3609.7.153",
        "id": "3609.7.153*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-8*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-8",
        "length": 312,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATGTCCTGTCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGGATATTGCAGCCCTCCCAGACCCTCAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTGAGCACTTTTGGTATGGGTGTAGGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTGGAGTGGCTGGCACACATTTGGTGGGATGATGATAAGTACTATAACCCAGCCCTGAAGAGTCGGCTCACAATCTCCAAGGATACCTCCAAAAACCAGGTATTCCTCAAGATCGCCAATGTGGACACTGCAGATACTGCCACATACTACTGTGCTCGAATAG",
            "read frame": 2,
            "strand": True
        },
        "start": 776704,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 666668,
                "length": 7,
                "sequence": {
                    "nucleotide": "TACATTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 666661
            },
            "3 nonamer": {
                "end": 666700,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACATTATTT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 666691
            },
            "end": 666661,
            "length": 303,
            "sequence": {
                "nucleotide": "GTTATTCTGAAAGAGTCTGGCCCTGGAATATTGCAGCCCTCCCAGACCCTAAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTTAGCACTTATGGTATGACAGTGAGCTGGATGCGTCAGCCTTCAAGGAAGGGTTTGGAGTGGCTGGCACAAATTGGGTCAGATGATAGTAAGCTCTATAACCCATCTCTGAAGAGCCGAATCACAATCTCCAAGGAAACCTCCAACAACCATGTATTCCTCAAGATCATCTGTGTGGACACTGAAGTTTCTGTCACATACTACTGTGCTCACAGACTACAG",
                "read frame": 2,
                "strand": True
            },
            "start": 666358
        },
        "accession": "BN000872",
        "allele": "3609.8pg.160*01",
        "confirmed": True,
        "end": 666661,
        "family": "3609v",
        "framed": True,
        "functionality": "O",
        "gene": "3609.8pg.160",
        "id": "3609.8pg.160*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-8-1*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-8-1",
        "length": 303,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTTATTCTGAAAGAGTCTGGCCCTGGAATATTGCAGCCCTCCCAGACCCTAAGTCTGACTTGTTCTTTCTCTGGGTTTTCACTTAGCACTTATGGTATGACAGTGAGCTGGATGCGTCAGCCTTCAAGGAAGGGTTTGGAGTGGCTGGCACAAATTGGGTCAGATGATAGTAAGCTCTATAACCCATCTCTGAAGAGCCGAATCACAATCTCCAAGGAAACCTCCAACAACCATGTATTCCTCAAGATCATCTGTGTGGACACTGAAGTTTCTGTCACATACTACTGTGCTCACAGACTACAG",
            "read frame": 2,
            "strand": True
        },
        "start": 666358,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 602746,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACGTTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 602739
            },
            "3 nonamer": {
                "end": 602778,
                "length": 9,
                "sequence": {
                    "nucleotide": "AAAGTATTT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 602769
            },
            "end": 602739,
            "length": 312,
            "sequence": {
                "nucleotide": "ATGTCCTATCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGTATATTGCAGCCCTCCCAGACCCTCAGTCTGACCTGTTCTTTCTCTGTGTTTTCACTGAGCACTTTTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTGGAGTGGCTGGCACACATTTATTGGGATGAGGACAAGCACTATAAACCATCCTTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCACTGTGGACACTGCAGATACTGCCACATACTACTCTGCTCGAAGAG",
                "read frame": 0,
                "strand": True
            },
            "start": 602427
        },
        "accession": "BN000872",
        "allele": "3609.9.164*01",
        "confirmed": True,
        "end": 602739,
        "family": "3609v",
        "framed": True,
        "functionality": "O",
        "gene": "3609.9.164",
        "id": "3609.9.164*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-9*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-9",
        "length": 312,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATGTCCTATCCCAGGTTACTCTGAAAGAGTCTGGCCCTGGTATATTGCAGCCCTCCCAGACCCTCAGTCTGACCTGTTCTTTCTCTGTGTTTTCACTGAGCACTTTTGGTATGGGTGTGAGCTGGATTCGTCAGCCTTCAGGGAAGGGTCTGGAGTGGCTGGCACACATTTATTGGGATGAGGACAAGCACTATAAACCATCCTTGAAGAGCCGGCTCACAATCTCCAAGGATACCTCCAACAACCAGGTATTCCTCAAGATCACCACTGTGGACACTGCAGATACTGCCACATACTACTCTGCTCGAAGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 602427,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 1738306,
            "length": 281,
            "sequence": {
                "nucleotide": "GTTTAGCTTCTAGAATCTGGAGGCTCTTTGATATAGCCTGGAGGATTCTTTAAACTCTCTTCTAAAGCTGCTGGATTCACCTTCAATGATTACTGGATGCACTGGTTTTGATAAGGCTCAGGGAATCGGCTAGAGTGGTTAGCAGATATAAAATATGTGAAAAGTGTAAAAAACCATGCAGAGTCTGTGAAGGGAAGATTTGCCACTACAGAGACAATTCTAAAAACTTATTTTATCTACAAATGAACAGTCTGAGAAATGAGGACATTGCCCCTTATTAC",
                "read frame": 1,
                "strand": True
            },
            "start": 1738025
        },
        "accession": "BN000872",
        "allele": "PG.15.75*01",
        "confirmed": True,
        "end": 1738313,
        "family": "3609v",
        "framed": False,
        "functionality": "P",
        "gene": "PG.15.75",
        "id": "PG.15.75*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV8-1*01",
        "imgt family name": "IGHV8",
        "imgt gene name": "IGHV8-1",
        "length": 291,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTTTAGCTTCTAGAATCTGGAGGCTCTTTGATATAGCCTGGAGGATTCTTTAAACTCTCTTCTAAAGCTGCTGGATTCACCTTCAATGATTACTGGATGCACTGGTTTTGATAAGGCTCAGGGAATCGGCTAGAGTGGTTAGCAGATATAAAATATGTGAAAAGTGTAAAAAACCATGCAGAGTCTGTGAAGGGAAGATTTGCCACTACAGAGACAATTCTAAAAACTTATTTTATCTACAAATGAACAGTCTGAGAAATGAGGACATTGCCCCTTATTACTCTATGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1738022,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2401575,
            "length": 299,
            "sequence": {
                "nucleotide": "GTGCCCAGTGGGAGGTAAAGCTGGTGGAGTCTAGGGGAGGCTTAGTGTAGCCTGGAAGGTCCGTGATACGCTCATGTGCAGCCTCTGGATTCACTGACAGTGACTAATAGTTGGCCTAGGTTTGCCAAGTTCCAAAGAAGGAGCTGGAATGGGGGGCATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGAAGAAGTGGGTTGGACATAACATATTCAGAATCAATATTTAGAGATTCTAATCCTTGAAGACATCACTTTTGACCAAGTATATATGCAACATGTTA",
                "read frame": 2,
                "strand": True
            },
            "start": 2401276
        },
        "accession": "BN000872",
        "allele": "7183.10pg.16*01",
        "confirmed": True,
        "end": 2401583,
        "family": "7183",
        "framed": False,
        "functionality": "P",
        "gene": "7183.10pg.16",
        "id": "7183.10pg.16*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-10*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-10",
        "length": 296,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTAAAGCTGGTGGAGTCTAGGGGAGGCTTAGTGTAGCCTGGAAGGTCCGTGATACGCTCATGTGCAGCCTCTGGATTCACTGACAGTGACTAATAGTTGGCCTAGGTTTGCCAAGTTCCAAAGAAGGAGCTGGAATGGGGGGCATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGAAGAAGTGGGTTGGACATAACATATTCAGAATCAATATTTAGAGATTCTAATCCTTGAAGACATCACTTTTGACCAAGTATATATGCAACATGTTACTTAGGTT",
            "read frame": 0,
            "strand": True
        },
        "start": 2401287,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2384029,
            "length": 278,
            "sequence": {
                "nucleotide": "GTGCCCAGTGTGAGGTGAAGCTGGTAGGGTCAGGGCAGCCTGGAGGGTCCCTGAAACACTCCTGTGCAGCCTCTGTAGTCACTGTGAGTGACTACTGAATGACCTGGGTCCTTCAGGCTCTAAAGAAGGGGCTGGAGAGGGTGGAAATAATTTTTAATGGTGGAGGTAGCACCTATTATCCAGACACCATGAAGGGCTGATTCACCATCTACAGAGATGATGCCAGAAACACACTTTACCTGAAAATAAACAGTCTGAGGTCTGAGTACACAGCCATG",
                "read frame": 1,
                "strand": True
            },
            "start": 2383751
        },
        "accession": "BN000872",
        "allele": "7183.11pg.19*01",
        "confirmed": True,
        "end": 2384027,
        "family": "7183",
        "framed": False,
        "functionality": "P",
        "gene": "7183.11pg.19",
        "id": "7183.11pg.19*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-11*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-11",
        "length": 265,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTAGGGTCAGGGCAGCCTGGAGGGTCCCTGAAACACTCCTGTGCAGCCTCTGTAGTCACTGTGAGTGACTACTGAATGACCTGGGTCCTTCAGGCTCTAAAGAAGGGGCTGGAGAGGGTGGAAATAATTTTTAATGGTGGAGGTAGCACCTATTATCCAGACACCATGAAGGGCTGATTCACCATCTACAGAGATGATGCCAGAAACACACTTTACCTGAAAATAAACAGTCTGAGGTCTGAGTACACAGCCA",
            "read frame": 0,
            "strand": True
        },
        "start": 2383762,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2369034,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAATG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2369027
            },
            "3 nonamer": {
                "end": 2369066,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACTAAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2369057
            },
            "end": 2369027,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGAAGTGAAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTGACTATTACATGTATTGGGTTCGCCAGACTCCAGAGAAGAGGCTGGAGTGGGTCGCATACATTAGTAATGGTGGTGGTAGCACCTATTATCCAGACACTGTAAAGGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCCGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGACA",
                "read frame": 2,
                "strand": True
            },
            "start": 2368720
        },
        "accession": "BN000872",
        "allele": "7183.12.20*01",
        "confirmed": True,
        "end": 2369027,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.12.20",
        "id": "7183.12.20*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-12*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-12",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTGAAGTGAAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTGACTATTACATGTATTGGGTTCGCCAGACTCCAGAGAAGAGGCTGGAGTGGGTCGCATACATTAGTAATGGTGGTGGTAGCACCTATTATCCAGACACTGTAAAGGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCCGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGACA",
            "read frame": 2,
            "strand": True
        },
        "start": 2368720,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2346396,
            "length": 296,
            "sequence": {
                "nucleotide": "GTGTTCAGGGTGAGGTAGAGCTGCTGAAATCCAAGGAAGGCATAGTGATACCTAAAAGCCCCTGAAACTCTCATGATAAACCTCTGGATTCACATTCAGTGATTACTATATAAGCAGGGTCCTCCACCCTCTAGGGAAGGGCCTTCAGTGGTGAGCATACATTAATAGAGACAGCTCCATCAACTGTGCTGAAGCCATGAAAAGCCAATTCACCATCACCAGAGACAATGCCAAGAACAACCTCTACCTGGATATTTGCAGACTGAGGACACAGCCATGTTTTCTTGTGCAGGATA",
                "read frame": 2,
                "strand": True
            },
            "start": 2346100
        },
        "accession": "BN000872",
        "allele": "7183.13pg.24*01",
        "confirmed": True,
        "end": 2346396,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.13pg.24",
        "id": "7183.13pg.24*01-C57BL/6",
        "identified": True,
        "imgt family name": "IGHV5",
        "length": 295,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGTTCAGGGTGAGGTAGAGCTGCTGAAATCCAAGGAAGGCATAGTGATACCTAAAAGCCCCTGAAACTCTCATGATAAACCTCTGGATTCACATTCAGTGATTACTATATAAGCAGGGTCCTCCACCCTCTAGGGAAGGGCCTTCAGTGGTGAGCATACATTAATAGAGACAGCTCCATCAACTGTGCTGAAGCCATGAAAAGCCAATTCACCATCACCAGAGACAATGCCAAGAACAACCTCTACCTGGATATTTGCAGACTGAGGACACAGCCATGTTTTCTTGTGCAGGATA",
            "read frame": 0,
            "strand": True
        },
        "start": 2346101,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2335048,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAATG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2335041
            },
            "3 nonamer": {
                "end": 2335080,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACTAAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2335071
            },
            "end": 2335041,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGACGTGAAGCTGGTGGAGTCTGGGGAAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGCCATGTCTTGGGTTCGCCAGACTCCAGAGAAGAGGCTGGAGTGGGTCGCATACATTAGTAGTGGTGGTGATTACATCTACTATGCAGACACTGTGAAGGGCCGATTCACCATCTCCAGAGACAATGCCAGGAACACCCTGTACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTACAAGAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 2334734
        },
        "accession": "BN000872",
        "allele": "7183.14.25*01",
        "confirmed": True,
        "end": 2335041,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.14.25",
        "id": "7183.14.25*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-9-1*02",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-9-1",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTGACGTGAAGCTGGTGGAGTCTGGGGAAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGCCATGTCTTGGGTTCGCCAGACTCCAGAGAAGAGGCTGGAGTGGGTCGCATACATTAGTAGTGGTGGTGATTACATCTACTATGCAGACACTGTGAAGGGCCGATTCACCATCTCCAGAGACAATGCCAGGAACACCCTGTACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTACAAGAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 2334734,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2325722,
            "length": 299,
            "sequence": {
                "nucleotide": "GTGCCCATTGTGAGGTGAAGCTGGTGGAGTCTAAGGGAGGTTTAGTGTAGCCTGGAAGGTCCATGATACTCTACTGTGCAGCCTCTGGATTCACTGTCAGTGACTACTGGTTGGCTGGGTTTGCCAGGCTCCAAAGAAGGAGCTGGAATGGGTGGGGACATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGTAAAGTGGGTTGGACATGACATATTCAGAATCAATATTTAGAGATTCTAATCCTTGAAGACATCACTTTCGACCAAGTATATATGCAAGATGTT",
                "read frame": 0,
                "strand": True
            },
            "start": 2325423
        },
        "accession": "BN000872",
        "allele": "7183.15pg.26*01",
        "confirmed": True,
        "end": 2325722,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.15pg.26",
        "id": "7183.15pg.26*01-C57BL/6",
        "identified": True,
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-7-6",
        "length": 298,
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGCCCATTGTGAGGTGAAGCTGGTGGAGTCTAAGGGAGGTTTAGTGTAGCCTGGAAGGTCCATGATACTCTACTGTGCAGCCTCTGGATTCACTGTCAGTGACTACTGGTTGGCTGGGTTTGCCAGGCTCCAAAGAAGGAGCTGGAATGGGTGGGGACATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGTAAAGTGGGTTGGACATGACATATTCAGAATCAATATTTAGAGATTCTAATCCTTGAAGACATCACTTTCGACCAAGTATATATGCAAGATGTT",
            "read frame": 0,
            "strand": True
        },
        "start": 2325424,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2308906,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2308899
            },
            "3 nonamer": {
                "end": 2308938,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2308929
            },
            "end": 2308899,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGAAGTGAAGCTCGTGGAGTCTGGGGGAGGTTTAGTGAAGCCTGGACAGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAACTATTACATGTCTTGGGTTCACCAGACTCCGGAGAAGAGGCTGGAGTGGGTTGCATACATTAGTAGTAGTGGTGTTAGCACCTATTATCCAGACAATGTAAAGGGCCGATTCGCCATCTCCAGAGACAATGCCAAGAACACCCTTTACCTGCAAATGACCAGTCTGAAGTCAGAGGACACGGCCTTGTATTACTGTGCAAGAGG",
                "read frame": 1,
                "strand": True
            },
            "start": 2308592
        },
        "accession": "BN000872",
        "allele": "7183.16.27*01",
        "confirmed": True,
        "end": 2308899,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.16.27",
        "id": "7183.16.27*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-12-4*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-12-4",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTGAAGTGAAGCTCGTGGAGTCTGGGGGAGGTTTAGTGAAGCCTGGACAGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAACTATTACATGTCTTGGGTTCACCAGACTCCGGAGAAGAGGCTGGAGTGGGTTGCATACATTAGTAGTAGTGGTGTTAGCACCTATTATCCAGACAATGTAAAGGGCCGATTCGCCATCTCCAGAGACAATGCCAAGAACACCCTTTACCTGCAAATGACCAGTCTGAAGTCAGAGGACACGGCCTTGTATTACTGTGCAAGAGG",
            "read frame": 1,
            "strand": True
        },
        "start": 2308592,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2289766,
            "length": 295,
            "sequence": {
                "nucleotide": "GTGTTCAGGGTGAGGTAGAGCTGCTGAAATCCAAGGAAGGCATAGTGATGCCTAAAAGTCCCTGAAACACTCATGTTAAACCTCTGGATTCACATTCAGTGATACTATATAAGCAGGGTCCTCCACCCTCTAGGGAAGGGTCTTCAGTTGTGAGCATACATTAATAGAGACAGCTCCATCAACTGTGCTGAAGCCATGAAAAGCCAACTCACCTTCACCAGAGACAATGCCAAGAACAACCTCTACCTGGAAATTTGCAGGCTGAGGACACAGCCATGTTTTCTTGTGCAGGATA",
                "read frame": 0,
                "strand": True
            },
            "start": 2289471
        },
        "accession": "BN000872",
        "allele": "7183.17pg.31*01",
        "confirmed": True,
        "end": 2289766,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.17pg.31",
        "id": "7183.17pg.31*01-C57BL/6",
        "identified": True,
        "imgt family name": "IGHV5",
        "length": 294,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGTTCAGGGTGAGGTAGAGCTGCTGAAATCCAAGGAAGGCATAGTGATGCCTAAAAGTCCCTGAAACACTCATGTTAAACCTCTGGATTCACATTCAGTGATACTATATAAGCAGGGTCCTCCACCCTCTAGGGAAGGGTCTTCAGTTGTGAGCATACATTAATAGAGACAGCTCCATCAACTGTGCTGAAGCCATGAAAAGCCAACTCACCTTCACCAGAGACAATGCCAAGAACAACCTCTACCTGGAAATTTGCAGGCTGAGGACACAGCCATGTTTTCTTGTGCAGGATA",
            "read frame": 0,
            "strand": True
        },
        "start": 2289472,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2244510,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2244503
            },
            "3 nonamer": {
                "end": 2244542,
                "length": 9,
                "sequence": {
                    "nucleotide": "GCACAAACT",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2244533
            },
            "end": 2244503,
            "length": 307,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGAGGTGAAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTGACTACGGAATGGCGTGGGTTCGACAGGCTCCAAGGAAGGGGCCTGAGTGGGTAGCATTCATTAGTAATTTGGCATATAGTATCTACTATGCAGACACTGTGACGGGCCGATTCACCATCTCTAGAGAGAATGCCAAGAACACCCTGTACCTGGAAATGAGCAGTCTGAGGTCTGAGGACACGGCCATGTATTACTGTGCAAGACA",
                "read frame": 2,
                "strand": True
            },
            "start": 2244196
        },
        "accession": "BN000872",
        "allele": "7183.18.35*01",
        "confirmed": True,
        "end": 2244503,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.18.35",
        "id": "7183.18.35*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-15*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-15",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTGAGGTGAAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTGACTACGGAATGGCGTGGGTTCGACAGGCTCCAAGGAAGGGGCCTGAGTGGGTAGCATTCATTAGTAATTTGGCATATAGTATCTACTATGCAGACACTGTGACGGGCCGATTCACCATCTCTAGAGAGAATGCCAAGAACACCCTGTACCTGGAAATGAGCAGTCTGAGGTCTGAGGACACGGCCATGTATTACTGTGCAAGACA",
            "read frame": 2,
            "strand": True
        },
        "start": 2244196,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2232630,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2232623
            },
            "3 nonamer": {
                "end": 2232662,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2232653
            },
            "end": 2232623,
            "length": 305,
            "sequence": {
                "nucleotide": "GAGTCCAGTGTGAAGTGAAGCTGGTGGAGTCTGAGGGAGGCTTAGTGCAGCCTGGAAGTTCCATGAAACTCTCCTGCACAGCCTCTGGATTCACTTTCAGTGACTATTACATGGCTTGGGTCCGCCAGGTTCCAGAAAAGGGTCTAGAATGGGTTGCAAACATTAATTATGATGGTAGTAGCACCTACTATCTGGACTCCTTGAAGAGCCGTTTCATCATCTCGAGAGACAATGCAAAGAACATTCTATACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCACGTATTACTGTGCAAGAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 2232316
        },
        "accession": "BN000872",
        "allele": "7183.19.36*01",
        "confirmed": True,
        "end": 2232623,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.19.36",
        "id": "7183.19.36*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-16*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-16",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGTCCAGTGTGAAGTGAAGCTGGTGGAGTCTGAGGGAGGCTTAGTGCAGCCTGGAAGTTCCATGAAACTCTCCTGCACAGCCTCTGGATTCACTTTCAGTGACTATTACATGGCTTGGGTCCGCCAGGTTCCAGAAAAGGGTCTAGAATGGGTTGCAAACATTAATTATGATGGTAGTAGCACCTACTATCTGGACTCCTTGAAGAGCCGTTTCATCATCTCGAGAGACAATGCAAAGAACATTCTATACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCACGTATTACTGTGCAAGAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 2232316,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "AF120460",
        "allele": "7183.1b*01",
        "confirmed": True,
        "end": 297,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.1b",
        "id": "7183.1b*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-2*03",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-2",
        "length": 297,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAAGTGAAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGAGTCCCTGAAACTCTCCTGTGAATCCAATGAATACGAATTCCCTTCCCATGACATGTCTTGGGTCCGCAAGACTCCGGAGAAGAGGCTGGAGTTGGTCGCAGCCATTAATAGTGATGGTGGTAGCACCTACTATCCAGACACCATGGAGAGACGATTCATCATCTCCAGAGACAATACCAAGAAGACCCTGTACCTGCAAATGAGCAGTCTGAGGTCTGAGGACACAGCCATGTATTACTGTGCAAGACGA",
            "read frame": 0,
            "strand": True
        },
        "start": 0,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2498336,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2498329
            },
            "3 nonamer": {
                "end": 2498368,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2498359
            },
            "end": 2498329,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGAAGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGGGTCCCGGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGCCATGTCTTGGGTCCGCCAGACTCCGGAGAAGAGGCTGGAGTGGGTCGCAGCCATTAGTACTGATGGTAGTTTCATCTACTAACCAGACACTGTAAAAGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTTTCTGCAAATGAGCAGTCTAAGGTATGAGGACACGGCCATGTATTACTGTTTGAGACA",
                "read frame": 0,
                "strand": True
            },
            "start": 2498022
        },
        "accession": "BN000872",
        "allele": "7183.1pg.1*01",
        "confirmed": True,
        "end": 2498329,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "7183.1pg.1",
        "id": "7183.1pg.1*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-1*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-1",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTGAAGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGGGTCCCGGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGCCATGTCTTGGGTCCGCCAGACTCCGGAGAAGAGGCTGGAGTGGGTCGCAGCCATTAGTACTGATGGTAGTTTCATCTACTAACCAGACACTGTAAAAGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTTTCTGCAAATGAGCAGTCTAAGGTATGAGGACACGGCCATGTATTACTGTTTGAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2498022,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2492759,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAATG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2492752
            },
            "3 nonamer": {
                "end": 2492791,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACTAAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2492782
            },
            "end": 2492752,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTACAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGAGTCCCTGAAACTCTCCTGTGAATCCAATGAATACGAATTCCCTTCCCATGACATGTCTTGGGTCCGCAAGACTCCGGAGAAGAGGCTGGAGTTGGTCGCAGCCATTAATAGTGATGGTGGTAGCACCTACTATCCAGACACCATGGAGAGACGATTCATCATCTCCAGAGACAATACCAAGAAGACCCTGTACCTGCAAATGAGCAGTCTGAGGTCTGAGGACACAGCCTTGTATTACTGTGCAAGACA",
                "read frame": 0,
                "strand": True
            },
            "start": 2492445
        },
        "accession": "BN000872",
        "allele": "7183.2.3*01",
        "confirmed": True,
        "end": 2492752,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.2.3",
        "id": "7183.2.3*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-2*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-2",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTACAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGCAGCCTGGAGAGTCCCTGAAACTCTCCTGTGAATCCAATGAATACGAATTCCCTTCCCATGACATGTCTTGGGTCCGCAAGACTCCGGAGAAGAGGCTGGAGTTGGTCGCAGCCATTAATAGTGATGGTGGTAGCACCTACTATCCAGACACCATGGAGAGACGATTCATCATCTCCAGAGACAATACCAAGAAGACCCTGTACCTGCAAATGAGCAGTCTGAGGTCTGAGGACACAGCCTTGTATTACTGTGCAAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2492445,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2212007,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2212000
            },
            "3 nonamer": {
                "end": 2212039,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2212030
            },
            "end": 2212000,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTGACTATGGAATGCACTGGGTTCGTCAGGCTCCAGAGAAGGGGCTGGAGTGGGTTGCATACATTAGTAGTGGCAGTAGTACCATCTACTATGCAGACACAGTGAAGGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTTCCTGCAAATGACCAGTCTGAGGTCTGAGGACACGGCCATGTATTACTGTGCAAGG",
                "read frame": 1,
                "strand": True
            },
            "start": 2211695
        },
        "accession": "BN000872",
        "allele": "7183.20.37*01",
        "confirmed": True,
        "end": 2212000,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.20.37",
        "id": "7183.20.37*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-17*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-17",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTGACTATGGAATGCACTGGGTTCGTCAGGCTCCAGAGAAGGGGCTGGAGTGGGTTGCATACATTAGTAGTGGCAGTAGTACCATCTACTATGCAGACACAGTGAAGGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTTCCTGCAAATGACCAGTCTGAGGTCTGAGGACACGGCCATGTATTACTGTGCAAGG",
            "read frame": 1,
            "strand": True
        },
        "start": 2211695,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2195377,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGCCCAGTGTGAGATGTAGCTGGTGGAGTCTAGGAGAGGCTTAGTACAGCCTGGAAGGTCCCTGAACTCTCATGAGCAGCCTCTGGATTCACTTTTAGTAATTATGTCATGGCCTGGGTCCAACAAGCTCCAGAGAAGGGGCTGGAGTGGGTTGCATACATTAGTAGTGGCAGTAGTACCTTCTATTATGCAGACACAGTTAAGGGCCCATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTTCCAGCAAATGAGCAGTCTAAGGTCTGAGGACAGAGACAAAAAGATAAAAATATAT",
                "read frame": 0,
                "strand": True
            },
            "start": 2195073
        },
        "accession": "BN000872",
        "allele": "7183.21pg.38*01",
        "confirmed": True,
        "end": 2195382,
        "family": "7183",
        "framed": False,
        "functionality": "P",
        "gene": "7183.21pg.38",
        "id": "7183.21pg.38*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-18*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-18",
        "length": 296,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GATGTAGCTGGTGGAGTCTAGGAGAGGCTTAGTACAGCCTGGAAGGTCCCTGAACTCTCATGAGCAGCCTCTGGATTCACTTTTAGTAATTATGTCATGGCCTGGGTCCAACAAGCTCCAGAGAAGGGGCTGGAGTGGGTTGCATACATTAGTAGTGGCAGTAGTACCTTCTATTATGCAGACACAGTTAAGGGCCCATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTTCCAGCAAATGAGCAGTCTAAGGTCTGAGGACAGAGACAAAAAGATAAAAATATATGGCCA",
            "read frame": 0,
            "strand": True
        },
        "start": 2195086,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2481291,
            "length": 278,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGAGGTGAAGCTGATAGGGTCAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGAGTCACTGTGAGTGACTACTGAATGACCTGGGTCCTTCAGGCTCCAAAGAAGGGGCTGGAGAGGGTGGCAATAATTTTTAATGGTGGAGGTAGCACCTATTATCCAGAAACCATGAAGGGCTGATTCACTATCTACAGAGATGATGCCACAAACACACTTTACCTGAAAATAAACAGTCTGAGGTCTTAGTACACAGCCCTG",
                "read frame": 2,
                "strand": True
            },
            "start": 2481013
        },
        "accession": "BN000872",
        "allele": "7183.3pg.5*01",
        "confirmed": True,
        "end": 2481289,
        "family": "7183",
        "framed": False,
        "functionality": "P",
        "gene": "7183.3pg.5",
        "id": "7183.3pg.5*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-3*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-3",
        "length": 265,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGATAGGGTCAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGAGTCACTGTGAGTGACTACTGAATGACCTGGGTCCTTCAGGCTCCAAAGAAGGGGCTGGAGAGGGTGGCAATAATTTTTAATGGTGGAGGTAGCACCTATTATCCAGAAACCATGAAGGGCTGATTCACTATCTACAGAGATGATGCCACAAACACACTTTACCTGAAAATAAACAGTCTGAGGTCTTAGTACACAGCCC",
            "read frame": 0,
            "strand": True
        },
        "start": 2481024,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2473817,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAATG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2473810
            },
            "3 nonamer": {
                "end": 2473849,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACTAAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2473840
            },
            "end": 2473810,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGAAGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGCCATGTCTTGGGTTCGCCAGACTCCGGAAAAGAGGCTGGAGTGGGTCGCAACCATTAGTGATGGTGGTAGTTACACCTACTATCCAGACAATGTAAAGGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACAACCTGTACCTGCAAATGAGCCATCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 2473503
        },
        "accession": "BN000872",
        "allele": "7183.4.6*01",
        "confirmed": True,
        "end": 2473810,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.4.6",
        "id": "7183.4.6*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-4*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-4",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTGAAGTGCAGCTGGTGGAGTCTGGGGGAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGCCATGTCTTGGGTTCGCCAGACTCCGGAAAAGAGGCTGGAGTGGGTCGCAACCATTAGTGATGGTGGTAGTTACACCTACTATCCAGACAATGTAAAGGGCCGATTCACCATCTCCAGAGACAATGCCAAGAACAACCTGTACCTGCAAATGAGCCATCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2473503,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2458262,
            "length": 294,
            "sequence": {
                "nucleotide": "CATTCTCTGCTTTGAGGTGTCAAGTGTGAGGTGAAGCTGGTGGAGTCTGTGGGAGGGTTAGAGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGAGTCACTGTCAGTGAGTACTGAATGACCTGGGTCCTTCAGGCTCCAAAGAAGGAGCTGGAAAGGGTGGCAATAATTTTTAATGGTGGAGGTAGCACCTATTATCCAGAAACCATGAAGGGCCAATTCACCATCTGCAGAGATGATGCCAAAAACATACTTTACCTAGAAATAAACACTCTGTGGTCTGAG",
                "read frame": 1,
                "strand": True
            },
            "start": 2457968
        },
        "accession": "BN000872",
        "allele": "7183.6pg.9*01",
        "confirmed": True,
        "end": 2458263,
        "family": "7183",
        "framed": False,
        "functionality": "P",
        "gene": "7183.6pg.9",
        "id": "7183.6pg.9*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-5*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-5",
        "length": 268,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTGGAGTCTGTGGGAGGGTTAGAGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGAGTCACTGTCAGTGAGTACTGAATGACCTGGGTCCTTCAGGCTCCAAAGAAGGAGCTGGAAAGGGTGGCAATAATTTTTAATGGTGGAGGTAGCACCTATTATCCAGAAACCATGAAGGGCCAATTCACCATCTGCAGAGATGATGCCAAAAACATACTTTACCTAGAAATAAACACTCTGTGGTCTGAGT",
            "read frame": 0,
            "strand": True
        },
        "start": 2457995,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2445757,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAATG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2445750
            },
            "3 nonamer": {
                "end": 2445789,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACTAAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2445780
            },
            "end": 2445750,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGACTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGGCATGTCTTGGGTTCGCCAGACTCCAGACAAGAGGCTGGAGTGGGTCGCAACCATTAGTAGTGGTGGTAGTTACACCTACTATCCAGACAGTGTGAAGGGGCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGACA",
                "read frame": 1,
                "strand": True
            },
            "start": 2445443
        },
        "accession": "BN000872",
        "allele": "7183.7.10*01",
        "confirmed": True,
        "end": 2445750,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.7.10",
        "id": "7183.7.10*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-6*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-6",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGACTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGGCATGTCTTGGGTTCGCCAGACTCCAGACAAGAGGCTGGAGTGGGTCGCAACCATTAGTAGTGGTGGTAGTTACACCTACTATCCAGACAGTGTGAAGGGGCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGACA",
            "read frame": 1,
            "strand": True
        },
        "start": 2445443,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2416289,
            "length": 288,
            "sequence": {
                "nucleotide": "GTGCCCAGTGTGAGGTGAAGCTGGTGGAGTCTGTGGGAGGCTTAGTGCAGCCTGGAGGGTCGCTGAAACTCTCCTGGGCAGCCTCTGGATTAACTCACTGACAACTGAATGACCTGGGTCCTTCAGGCTCCAAGGAAGGGGCTGGAGAGGGTGGCAATAATTTTTAGTGGTGGAGGTAGCACCTACTATATAGACGCAATGAAGGGCCGATTCACCGTCTCCAGAGATGATAACAAAAACACACTTTACCTGAAAATAAACAGTCTGAGGTCTGAGTACACAGCCATG",
                "read frame": 1,
                "strand": True
            },
            "start": 2416001
        },
        "accession": "BN000872",
        "allele": "7183.8pg.14*01",
        "confirmed": True,
        "end": 2416289,
        "family": "7183",
        "framed": False,
        "functionality": "P",
        "gene": "7183.8pg.14",
        "id": "7183.8pg.14*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-8*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-8",
        "length": 277,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTGGAGTCTGTGGGAGGCTTAGTGCAGCCTGGAGGGTCGCTGAAACTCTCCTGGGCAGCCTCTGGATTAACTCACTGACAACTGAATGACCTGGGTCCTTCAGGCTCCAAGGAAGGGGCTGGAGAGGGTGGCAATAATTTTTAGTGGTGGAGGTAGCACCTACTATATAGACGCAATGAAGGGCCGATTCACCGTCTCCAGAGATGATAACAAAAACACACTTTACCTGAAAATAAACAGTCTGAGGTCTGAGTACACAGCCATG",
            "read frame": 0,
            "strand": True
        },
        "start": 2416012,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2409494,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2409487
            },
            "3 nonamer": {
                "end": 2409526,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACTAAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2409517
            },
            "end": 2409487,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGAAGTGATGCTGGTGGAGTCTGGGGGAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATACCATGTCTTGGGTTCGCCAGACTCCGGAGAAGAGGCTGGAGTGGGTCGCAACCATTAGTGGTGGTGGTGGTAACACCTACTATCCAGACAGTGTGAAGGGTCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCAGTCTGAGGTCTGAGGACACGGCCTTGTATTACTGTGCAAGACA",
                "read frame": 0,
                "strand": True
            },
            "start": 2409180
        },
        "accession": "BN000872",
        "allele": "7183.9.15*01",
        "confirmed": True,
        "end": 2409487,
        "family": "7183",
        "framed": True,
        "functionality": "F",
        "gene": "7183.9.15",
        "id": "7183.9.15*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-9*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-9",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTGAAGTGATGCTGGTGGAGTCTGGGGGAGGCTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATACCATGTCTTGGGTTCGCCAGACTCCGGAGAAGAGGCTGGAGTGGGTCGCAACCATTAGTGGTGGTGGTGGTAACACCTACTATCCAGACAGTGTGAAGGGTCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCAGTCTGAGGTCTGAGGACACGGCCTTGTATTACTGTGCAAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2409180,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2436563,
            "length": 300,
            "sequence": {
                "nucleotide": "GTGCCCAGTGGGAGGTGAAGCTGGTGGAGTCTAGGGGAGGCTTAGTGTAGCCTGGAAGGTCCGTGATATGCTCATGTGCAGCCTCTGGATTCACTGACAGTGACTAATAGTTGGCCTAGGTTTGCCAAGTTCCAAAGAAGGAGCTGGAATGGGGGGCATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGAAGAAGTGGGTTGGACATAACATATTCAGAATCAATATTTAGAGATTCTAATCCTTGAAGACATCACTTTTGACCAAGTATATATGCAACATGTTAC",
                "read frame": 1,
                "strand": True
            },
            "start": 2436263
        },
        "accession": "BN000872",
        "allele": "PG.1.11*01",
        "confirmed": True,
        "end": 2436570,
        "family": "7183",
        "framed": False,
        "functionality": "P",
        "gene": "PG.1.11",
        "id": "PG.1.11*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-7*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-7",
        "length": 296,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTGGTGGAGTCTAGGGGAGGCTTAGTGTAGCCTGGAAGGTCCGTGATATGCTCATGTGCAGCCTCTGGATTCACTGACAGTGACTAATAGTTGGCCTAGGTTTGCCAAGTTCCAAAGAAGGAGCTGGAATGGGGGGCATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGAAGAAGTGGGTTGGACATAACATATTCAGAATCAATATTTAGAGATTCTAATCCTTGAAGACATCACTTTTGACCAAGTATATATGCAACATGTTACTTAGGTT",
            "read frame": 0,
            "strand": True
        },
        "start": 2436274,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1751087,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1751080
            },
            "3 nonamer": {
                "end": 1751119,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAAAC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1751110
            },
            "end": 1751080,
            "length": 297,
            "sequence": {
                "nucleotide": "GGTGAGGTGAGGCTGGTGGAATCTGGAGGCAGCTTGATAGAGCATGAAGGGTACATTCAACTCTTTTGTCAAGCTTCTGGATTCACCCTCAGTGGTTACTGGATGCACTGGATTTGCCAAGCCCCAGGGAAGGGGCCAGAGTGGTTAGCAAATATAAAATATGATGGGAGTGAAAAATACTATGCAGTGTCTATGAAGGGGTGATTTGCCATCTCCAGAGACCTTCCTAAGAACTTTCTTTATCTGCAAATGAGCAATTTGAGAAATGAGGACACTGCAATGTATTACTGTGCAAGAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1750781
        },
        "accession": "BN000872",
        "allele": "PG.14.73*01",
        "confirmed": True,
        "end": 1751080,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "PG.14.73",
        "id": "PG.14.73*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-21*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-21",
        "length": 299,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGTGAGGTGAGGCTGGTGGAATCTGGAGGCAGCTTGATAGAGCATGAAGGGTACATTCAACTCTTTTGTCAAGCTTCTGGATTCACCCTCAGTGGTTACTGGATGCACTGGATTTGCCAAGCCCCAGGGAAGGGGCCAGAGTGGTTAGCAAATATAAAATATGATGGGAGTGAAAAATACTATGCAGTGTCTATGAAGGGGTGATTTGCCATCTCCAGAGACCTTCCTAAGAACTTTCTTTATCTGCAAATGAGCAATTTGAGAAATGAGGACACTGCAATGTATTACTGTGCAAGAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1750781,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2356544,
            "length": 290,
            "sequence": {
                "nucleotide": "TGTGAGGTGAAGCTGGTGAAGTCTAAGGGGAGGCATAGTGCAGCCTAGAAGGTCCATGATACTCTACTGTGCAGCCTCGGATTCACTGTAAGTGACGACTGGTTTGTCCGTGTTTGCCAGGCTCCAAAGAAGGGGCTGCAGTGGGGGATGGGAATAATTTTTCATGGTTGTGGTAGCCCCTCTTATGCAGACACCTTGAAGAAGTGGGTTGGACGTAACATATTCAGAATCAATATTTAAAGATTCTAATCCTTGAAGATATCACTTTTGACCAAGTATATATGAACCAT",
                "read frame": 0,
                "strand": True
            },
            "start": 2356254
        },
        "accession": "BN000872",
        "allele": "PG.2.21*01",
        "confirmed": True,
        "end": 2356544,
        "family": "7183",
        "framed": True,
        "functionality": "P",
        "gene": "PG.2.21",
        "id": "PG.2.21*01-C57BL/6",
        "identified": True,
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-13",
        "length": 289,
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGAGGTGAAGCTGGTGAAGTCTAAGGGGAGGCATAGTGCAGCCTAGAAGGTCCATGATACTCTACTGTGCAGCCTCGGATTCACTGTAAGTGACGACTGGTTTGTCCGTGTTTGCCAGGCTCCAAAGAAGGGGCTGCAGTGGGGGATGGGAATAATTTTTCATGGTTGTGGTAGCCCCTCTTATGCAGACACCTTGAAGAAGTGGGTTGGACGTAACATATTCAGAATCAATATTTAAAGATTCTAATCCTTGAAGATATCACTTTTGACCAAGTATATATGAACCAT",
            "read frame": 0,
            "strand": True
        },
        "start": 2356255,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "end": 2186386,
            "length": 297,
            "sequence": {
                "nucleotide": "ATATCCAGTGTGAGGTGCAGCTGGTGGGGTCCAGGAGAGTCTTAATGCAGCCTAGAAAAACCTTGAAACTCCCCTCTGCAGCACATGGAATCACCTTCAGCGACTGACATGGTACTCCGGGCTCCAGGGAAGGGGCTGGAGTGGGATACATTCATTAATAGTAAGGGTAGTTACATCAACTATGCCAAAGCTATAAAAGAACAATTCACTATCTCCAGAGACAATACCAAGAACACTTTGTACCTGGAAATTAGCAGTTTGAAGTCTGAAAACACAGTCATTTATTACTGTACAACA",
                "read frame": 2,
                "strand": True
            },
            "start": 2186089
        },
        "accession": "BN000872",
        "allele": "PG.7.41*01",
        "confirmed": True,
        "end": 2186386,
        "family": "7183",
        "framed": False,
        "functionality": "P",
        "gene": "PG.7.41",
        "id": "PG.7.41*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV5-19*01",
        "imgt family name": "IGHV5",
        "imgt gene name": "IGHV5-19",
        "length": 297,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATATCCAGTGTGAGGTGCAGCTGGTGGGGTCCAGGAGAGTCTTAATGCAGCCTAGAAAAACCTTGAAACTCCCCTCTGCAGCACATGGAATCACCTTCAGCGACTGACATGGTACTCCGGGCTCCAGGGAAGGGGCTGGAGTGGGATACATTCATTAATAGTAAGGGTAGTTACATCAACTATGCCAAAGCTATAAAAGAACAATTCACTATCTCCAGAGACAATACCAAGAACACTTTGTACCTGGAAATTAGCAGTTTGAAGTCTGAAAACACAGTCATTTATTACTGTACAACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2186089,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1602041,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1602034
            },
            "3 nonamer": {
                "end": 1602073,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1602064
            },
            "end": 1602034,
            "length": 332,
            "sequence": {
                "nucleotide": "GTGTCTACTCACAGGTTCAGCTCCAGCAGTCTGGGCCTGAGCTGGCAAGGCCTTGGGCTTCAGTGAAGATATCCTGCCAGGCTTTCTACACCTTTTCCAGAAGGGTGCACTTTGCCATTAGGGATACCAACTACTGGATGCAGTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATCGGGGCTATTTATCCTGGAAATGGTGATACTAGTTACAATCAGAAGTTCAAGGGCAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAACTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCATGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1601702
        },
        "accession": "BN000872",
        "allele": "J558.1.85*01",
        "confirmed": True,
        "end": 1602034,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.1.85",
        "id": "J558.1.85*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-2*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-2",
        "length": 332,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCTACTCACAGGTTCAGCTCCAGCAGTCTGGGCCTGAGCTGGCAAGGCCTTGGGCTTCAGTGAAGATATCCTGCCAGGCTTTCTACACCTTTTCCAGAAGGGTGCACTTTGCCATTAGGGATACCAACTACTGGATGCAGTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATCGGGGCTATTTATCCTGGAAATGGTGATACTAGTTACAATCAGAAGTTCAAGGGCAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAACTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCATGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1601702,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1440435,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1440428
            },
            "3 nonamer": {
                "end": 1440467,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1440458
            },
            "end": 1440428,
            "length": 291,
            "sequence": {
                "nucleotide": "GTCCAACTTCAGCAGTCTGGACCAGAGCTGGTAATACCTGGGGCTTAGGTGAAGTTGTCCTGCAAGGCTTCTGGCTACAATTTTAATGACTATGAAATTCAATGGGTGAAGCAGAGTCTGAAGCAGGGACTGGAATGGATTGGAGCTATTCATCCTGAAAATGGTGGTATTACCTACAATCAGAAGTTCAAAGGCAAGGCCACATTTACTGTAGACACATCCTCCAACACAGCCTACATGCAACTCAGAAGCCTGACATCTGAGGACACTGCTGACTATTATTGTGAGAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1440137
        },
        "accession": "BN000872",
        "allele": "J558.10pg.100*01",
        "confirmed": True,
        "end": 1440428,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.10pg.100",
        "id": "J558.10pg.100*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-13*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-13",
        "length": 291,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTCCAACTTCAGCAGTCTGGACCAGAGCTGGTAATACCTGGGGCTTAGGTGAAGTTGTCCTGCAAGGCTTCTGGCTACAATTTTAATGACTATGAAATTCAATGGGTGAAGCAGAGTCTGAAGCAGGGACTGGAATGGATTGGAGCTATTCATCCTGAAAATGGTGGTATTACCTACAATCAGAAGTTCAAAGGCAAGGCCACATTTACTGTAGACACATCCTCCAACACAGCCTACATGCAACTCAGAAGCCTGACATCTGAGGACACTGCTGACTATTATTGTGAGAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1440137,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1424543,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1424536
            },
            "3 nonamer": {
                "end": 1424575,
                "length": 9,
                "sequence": {
                    "nucleotide": "CATAAACCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1424566
            },
            "end": 1424536,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGTTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACATTCACTAGCTATGTTATGCACTGGGTGAAGCAGAAGCCTGGGCAGGGCCTTGAGTGGATTGGATATATTTATCCTTACAATGATGGTACTAAGTACAATGAGAAGTTCAAAGGCAAGGCCACACTGACTTCAGACAAATCCTCCAGCACAGCCTACATGGAGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1424231
        },
        "accession": "BN000872",
        "allele": "J558.11pg.101*01",
        "confirmed": True,
        "end": 1424536,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.11pg.101",
        "id": "J558.11pg.101*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-14*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-14",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGTTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACATTCACTAGCTATGTTATGCACTGGGTGAAGCAGAAGCCTGGGCAGGGCCTTGAGTGGATTGGATATATTTATCCTTACAATGATGGTACTAAGTACAATGAGAAGTTCAAAGGCAAGGCCACACTGACTTCAGACAAATCCTCCAGCACAGCCTACATGGAGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1424231,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1413762,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1413755
            },
            "3 nonamer": {
                "end": 1413794,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1413785
            },
            "end": 1413755,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAATCCCAGGTTCAACTGCAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGACGCTGTCCTGCAAGGCTTCGGGCTACACATTTACTGACTATGAAATGCACTGGGTGAAGCAGACACCTGTGCATGGCCTGGAATGGATTGGAGCTATTGATCCTGAAACTGGTGGTACTGCCTACAATCAGAAGTTCAAGGGCAAGGCCATACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCCGTCTATTACTGTACAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1413450
        },
        "accession": "BN000872",
        "allele": "J558.12.102*01",
        "confirmed": True,
        "end": 1413755,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.12.102",
        "id": "J558.12.102*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-15*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-15",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAATCCCAGGTTCAACTGCAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGACGCTGTCCTGCAAGGCTTCGGGCTACACATTTACTGACTATGAAATGCACTGGGTGAAGCAGACACCTGTGCATGGCCTGGAATGGATTGGAGCTATTGATCCTGAAACTGGTGGTACTGCCTACAATCAGAAGTTCAAGGGCAAGGCCATACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCCGTCTATTACTGTACAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1413450,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1405243,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1405236
            },
            "3 nonamer": {
                "end": 1405274,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1405265
            },
            "end": 1405236,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTAAGGTAGTGAATGCTGGGGCTTCCGTGAAGCTGTCCTGCAAGTCTTCTGGTTACTCATTCAGTAGATACAAAATGGAATGTGTGAAACAGAGCCATGTAAAGAGCCTTGAGTGGATTGAACATATTAATCTTTTCAATGGTATTACTAACTACAATGGAAACTTTAAAAGCAAGGCCACATTGACTGTAGACATATCCTCTAGCACAGCCTATATGGAGCTTAGCAGATTGACATCTGAAGACTCAGAGGTATATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1404931
        },
        "accession": "BN000872",
        "allele": "J558.13.103*01",
        "confirmed": True,
        "end": 1405236,
        "family": "J558",
        "framed": True,
        "functionality": "O",
        "gene": "J558.13.103",
        "id": "J558.13.103*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-16*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-16",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTAAGGTAGTGAATGCTGGGGCTTCCGTGAAGCTGTCCTGCAAGTCTTCTGGTTACTCATTCAGTAGATACAAAATGGAATGTGTGAAACAGAGCCATGTAAAGAGCCTTGAGTGGATTGAACATATTAATCTTTTCAATGGTATTACTAACTACAATGGAAACTTTAAAAGCAAGGCCACATTGACTGTAGACATATCCTCTAGCACAGCCTATATGGAGCTTAGCAGATTGACATCTGAAGACTCAGAGGTATATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1404931,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1399147,
                "length": 7,
                "sequence": {
                    "nucleotide": "TACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1399140
            },
            "3 nonamer": {
                "end": 1399178,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAATC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1399169
            },
            "end": 1399140,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCAACTCTGAAGTCCACCTGCAGCAGTCTCTACCTAAGGTAGTGAAGGCTGGGCCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGATCCTTCAAAGGATTGAATATGTTAATCCTTATAATGGTGGTACTGGCTACAATGAAAAGTTCAAGGACAAGGCCACATTGACTGCAGACAAATCCTTCAGCACAGCCTATATGCAATTCAGCAGCCTGACATCTGAGGACTCTCTGGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1398835
        },
        "accession": "BN000872",
        "allele": "J558.14pg.104*01",
        "confirmed": True,
        "end": 1399140,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.14pg.104",
        "id": "J558.14pg.104*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-17-1*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-17-1",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCAACTCTGAAGTCCACCTGCAGCAGTCTCTACCTAAGGTAGTGAAGGCTGGGCCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGATCCTTCAAAGGATTGAATATGTTAATCCTTATAATGGTGGTACTGGCTACAATGAAAAGTTCAAGGACAAGGCCACATTGACTGCAGACAAATCCTTCAGCACAGCCTATATGCAATTCAGCAGCCTGACATCTGAGGACTCTCTGGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1398835,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1394366,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1394359
            },
            "3 nonamer": {
                "end": 1394398,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1394389
            },
            "end": 1394359,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTATTCATTCAGTGACTACTACATGGAATGGGTGAAGCAGAGCCATAGAAAGAGCCTTGAATGTATTGGAGAAATTAATCCTTACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCTTCCAGCACAGCGTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTTCGGTCTATTACTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1394055
        },
        "accession": "BN000872",
        "allele": "J558.15pg.105*01",
        "confirmed": True,
        "end": 1394359,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.15pg.105",
        "id": "J558.15pg.105*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-17*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-17",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTATTCATTCAGTGACTACTACATGGAATGGGTGAAGCAGAGCCATAGAAAGAGCCTTGAATGTATTGGAGAAATTAATCCTTACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCTTCCAGCACAGCGTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTTCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1394055,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1388483,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1388476
            },
            "3 nonamer": {
                "end": 1388515,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1388506
            },
            "end": 1388476,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATACCCTGCAAGGCTTCTGGATACACATTCACTGACTACAACATGGACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGATATTAATCCTAACAATGGTGGTACTATCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACACTGCAGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1388171
        },
        "accession": "BN000872",
        "allele": "J558.16.106*01",
        "confirmed": True,
        "end": 1388476,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.16.106",
        "id": "J558.16.106*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-18*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-18",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATACCCTGCAAGGCTTCTGGATACACATTCACTGACTACAACATGGACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGATATTAATCCTAACAATGGTGGTACTATCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACACTGCAGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1388171,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1366940,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1366933
            },
            "3 nonamer": {
                "end": 1366972,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAGAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1366963
            },
            "end": 1366933,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGAAGCAGCCTGGAACTGTGGTGGTGAAACCTGGGGCTTCAGTGAAGATATCCTGCCAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTGTAGTGGATTGGACTTATTATTCCTTACAATGGTGATACTGGCTACAACCAGAAGTTCAAGGGCAAAGCAACATTGACTGTAGACAAGTCCTCCAGCACAGCCAACATGGAGCTCTGCAGCCTGACATCTGAGGACTCTACAGTCTATTACCGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1366628
        },
        "accession": "BN000872",
        "allele": "J558.17pg.107*01",
        "confirmed": True,
        "end": 1366933,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.17pg.107",
        "id": "J558.17pg.107*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-19-1*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-19-1",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGAAGCAGCCTGGAACTGTGGTGGTGAAACCTGGGGCTTCAGTGAAGATATCCTGCCAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTGTAGTGGATTGGACTTATTATTCCTTACAATGGTGATACTGGCTACAACCAGAAGTTCAAGGGCAAAGCAACATTGACTGTAGACAAGTCCTCCAGCACAGCCAACATGGAGCTCTGCAGCCTGACATCTGAGGACTCTACAGTCTATTACCGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1366628,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1362439,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1362432
            },
            "3 nonamer": {
                "end": 1362470,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAAAC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1362461
            },
            "end": 1362432,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAACAGTCTGGACCTGTGCTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGTAAGGCTTCTGGATACACATTCACTGACTACTATATGAACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGTTATTAATCCTTACAACGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTTGACAAGTCCTCCAGCACAGCCTACATGGAGCTCAACAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1362127
        },
        "accession": "BN000872",
        "allele": "J558.18.108*01",
        "confirmed": True,
        "end": 1362432,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.18.108",
        "id": "J558.18.108*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-19*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-19",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAACAGTCTGGACCTGTGCTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGTAAGGCTTCTGGATACACATTCACTGACTACTATATGAACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGTTATTAATCCTTACAACGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTTGACAAGTCCTCCAGCACAGCCTACATGGAGCTCAACAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1362127,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1347315,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACTGTA",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1347308
            },
            "3 nonamer": {
                "end": 1347347,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACT",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1347338
            },
            "end": 1347308,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTGTTCTCTGAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGATTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTTACTGGCTACTTTATGAACTGGGTGATGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACGTATTAATCCTTACAATGGTGATACTTTCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCTAGCACAGCCCACATGGAGCTCCGGAGCCTGACATCTGAGGACTCTGCAGTCTATTATTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1347003
        },
        "accession": "BN000872",
        "allele": "J558.19.109*01",
        "confirmed": True,
        "end": 1347308,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.19.109",
        "id": "J558.19.109*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-20*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-20",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTGTTCTCTGAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGATTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTTACTGGCTACTTTATGAACTGGGTGATGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACGTATTAATCCTTACAATGGTGATACTTTCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCTAGCACAGCCCACATGGAGCTCCGGAGCCTGACATCTGAGGACTCTGCAGTCTATTATTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1347003,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1583980,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1583973
            },
            "3 nonamer": {
                "end": 1584011,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1584002
            },
            "end": 1583973,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGGGCTGAACTGGCAAGACCTGGTGCCTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTTACTAGCTACACGATGCACTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATTGGATACATTAATCCTAGCAGTGGTTATACTAAGTACAATCAGAAGTTCAAGGACAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAACTGAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1583668
        },
        "accession": "BN000872",
        "allele": "J558.2.88*01",
        "confirmed": True,
        "end": 1583973,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.2.88",
        "id": "J558.2.88*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-4*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-4",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGGGCTGAACTGGCAAGACCTGGTGCCTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTTACTAGCTACACGATGCACTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATTGGATACATTAATCCTAGCAGTGGTTATACTAAGTACAATCAGAAGTTCAAGGACAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAACTGAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1583668,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1334842,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1334835
            },
            "3 nonamer": {
                "end": 1334873,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1334864
            },
            "end": 1334835,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAACAGTCTGGACCTGAGTTGGTGAAGCCTGGGGCTTTAGTGAAGTTATCCTGCAAGGTTTCTGGATTCACATTCACTGACTAATACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACATGTTTATCCTTACAATGGTGGTACTAGCTACAACCAGAAATTCAAGGGCAAGGCCACATTGACTGTCGACAATACCTCCAGCACAGCCTACATGGAGCTCGGCAGCCTGACTTCTGAGGACTCTGCGGTCTATTACTCTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1334530
        },
        "accession": "BN000872",
        "allele": "J558.20pg.110*01",
        "confirmed": True,
        "end": 1334835,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.20pg.110",
        "id": "J558.20pg.110*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-21-1*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-21-1",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAACAGTCTGGACCTGAGTTGGTGAAGCCTGGGGCTTTAGTGAAGTTATCCTGCAAGGTTTCTGGATTCACATTCACTGACTAATACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACATGTTTATCCTTACAATGGTGGTACTAGCTACAACCAGAAATTCAAGGGCAAGGCCACATTGACTGTCGACAATACCTCCAGCACAGCCTACATGGAGCTCGGCAGCCTGACTTCTGAGGACTCTGCGGTCTATTACTCTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1334530,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1330679,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1330672
            },
            "3 nonamer": {
                "end": 1330711,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1330702
            },
            "end": 1330672,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCATGGCTTCTGGTTATTCATTCAGTGACTACTACATGCACTGAGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGATATATTAACCCTAACAATGGTTGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCGTACATGGAGCTCCACAGCCTGACATCTGAGGACTCTTCGGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1330367
        },
        "accession": "BN000872",
        "allele": "J558.21pg.111*01",
        "confirmed": True,
        "end": 1330672,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.21pg.111",
        "id": "J558.21pg.111*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-21*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-21",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCATGGCTTCTGGTTATTCATTCAGTGACTACTACATGCACTGAGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGATATATTAACCCTAACAATGGTTGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCGTACATGGAGCTCCACAGCCTGACATCTGAGGACTCTTCGGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1330367,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1324814,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1324807
            },
            "3 nonamer": {
                "end": 1324846,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAATCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1324837
            },
            "end": 1324807,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACATTCACTGACTACAACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGATATATTAACCCTAACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAAACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCGGAGGATTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1324502
        },
        "accession": "BN000872",
        "allele": "J558.22.112*01",
        "confirmed": True,
        "end": 1324807,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.22.112",
        "id": "J558.22.112*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-22*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-22",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACATTCACTGACTACAACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGATATATTAACCCTAACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAAACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCGGAGGATTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1324502,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1306637,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1306630
            },
            "3 nonamer": {
                "end": 1306669,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1306660
            },
            "end": 1306630,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAATCCCAGGTTCAACTGCAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCGGGCTACACATTTACTGACTATGAAATGCACTGTGTGAAGCAGACACCTGTGCACGGCCTGGAATGGATTGGAGCTATTGATCCTGAAACTTGTGGTACTGCCTACAATCAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCCGTCTATTACTGTACAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1306325
        },
        "accession": "BN000872",
        "allele": "J558.23.113*01",
        "confirmed": True,
        "end": 1306630,
        "family": "J558",
        "framed": True,
        "functionality": "O",
        "gene": "J558.23.113",
        "id": "J558.23.113*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-23*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-23",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAATCCCAGGTTCAACTGCAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCGGGCTACACATTTACTGACTATGAAATGCACTGTGTGAAGCAGACACCTGTGCACGGCCTGGAATGGATTGGAGCTATTGATCCTGAAACTTGTGGTACTGCCTACAATCAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCCGTCTATTACTGTACAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1306325,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 nonamer": {
                "end": 1298188,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1298179
            },
            "end": 1298152,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCGGCTGCAGCAGTCTGGACCTAAGGTAGTGAATGCTGGGGCTTCCGTGAAGCTGTCCTGCAAGTCTTCTGGTTACTCATTCAGTAGATACAAAATGGAATGTGTGAAACAGAGCCATGGAAAGAGCCTTGAGTGGATTGAACATATTAATCTTTTCAAAGGTGTTACTAACTACAATGGAAAGTTCAAAAGCAAGGCCACATTGACTGTAGACATATCCTCTAGCACAGCCTATATGGAGCTTAGCAGATTGACATCTGAAGACTCAGAGGTATATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1297847
        },
        "accession": "BN000872",
        "allele": "J558.24pg.114*01",
        "confirmed": True,
        "end": 1298152,
        "family": "J558",
        "framed": True,
        "functionality": "O",
        "gene": "J558.24pg.114",
        "id": "J558.24pg.114*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-24*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-24",
        "length": 294,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTCCGGCTGCAGCAGTCTGGACCTAAGGTAGTGAATGCTGGGGCTTCCGTGAAGCTGTCCTGCAAGTCTTCTGGTTACTCATTCAGTAGATACAAAATGGAATGTGTGAAACAGAGCCATGGAAAGAGCCTTGAGTGGATTGAACATATTAATCTTTTCAAAGGTGTTACTAACTACAATGGAAAGTTCAAAAGCAAGGCCACATTGACTGTAGACATATCCTCTAGCACAGCCTATATGGAGCTTAGCAGATTGACATCTGAAGACTCAGAGGTATATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1297858,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1289678,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1289671
            },
            "3 nonamer": {
                "end": 1289710,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1289701
            },
            "end": 1289671,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCGTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACATTATGAACTTGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTTGAGAAATTAATCCTTACAATGGTGGTACTAACTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGAGTACATGGAGCTCCACAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1289366
        },
        "accession": "BN000872",
        "allele": "J558.25pg.115*01",
        "confirmed": True,
        "end": 1289671,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.25pg.115",
        "id": "J558.25pg.115*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-25*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-25",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCGTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACATTATGAACTTGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTTGAGAAATTAATCCTTACAATGGTGGTACTAACTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGAGTACATGGAGCTCCACAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1289366,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1282715,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1282708
            },
            "3 nonamer": {
                "end": 1282747,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1282738
            },
            "end": 1282708,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAATCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGTAAGGCTTCTGGATACACGTTCACTGACTACTACATGAACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGATATTAATCCTAACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1282403
        },
        "accession": "BN000872",
        "allele": "J558.26.116*01",
        "confirmed": True,
        "end": 1282708,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.26.116",
        "id": "J558.26.116*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-26*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-26",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAATCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGTAAGGCTTCTGGATACACGTTCACTGACTACTACATGAACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGATATTAATCCTAACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1282403,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1266471,
                "length": 7,
                "sequence": {
                    "nucleotide": "CAAAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1266464
            },
            "3 nonamer": {
                "end": 1266503,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1266494
            },
            "end": 1266464,
            "length": 305,
            "sequence": {
                "nucleotide": "GCATCCACTCTAAGGTCAAGCTGTAGCAGTCTGGACCTGAGCTGGTGAAGTCTGGGGCTTCAGAGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGATCATTGAGTAGATTGGACTTATTATTCCCTACAATGGTGATACTGGCTACAACCAGAAGTTTAAGGGCAAGGCCACATTGACTGTAGACTAGTCCTTCAGCACAGCCTACATGGAGCTCCACAGCCTGAAATCTTAGGACTCTGTGGTCTATTTCTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1266159
        },
        "accession": "BN000872",
        "allele": "J558.27pg.117*01",
        "confirmed": True,
        "end": 1266464,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.27pg.117",
        "id": "J558.27pg.117*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-27*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-27",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GCATCCACTCTAAGGTCAAGCTGTAGCAGTCTGGACCTGAGCTGGTGAAGTCTGGGGCTTCAGAGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGATCATTGAGTAGATTGGACTTATTATTCCCTACAATGGTGATACTGGCTACAACCAGAAGTTTAAGGGCAAGGCCACATTGACTGTAGACTAGTCCTTCAGCACAGCCTACATGGAGCTCCACAGCCTGAAATCTTAGGACTCTGTGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1266159,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1261666,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACGGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1261659
            },
            "3 nonamer": {
                "end": 1261697,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAAAC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1261688
            },
            "end": 1261659,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGATCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGTGCTTCAGTGAAGATATCCTTCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGAACTGGATGAAGTAGAGCCACGGAAATAGCCTTGAGTGGATTGGATATATTGATCCTTACAATGGTGGTACTAGCTACAACCAGAAATTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCTTCCAGCACAGCCTACATGCATCTCAACAACCTGACATCTGAGGACTCTGCAGTCTATTACTTTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1261354
        },
        "accession": "BN000872",
        "allele": "J558.28pg.118*01",
        "confirmed": True,
        "end": 1261659,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.28pg.118",
        "id": "J558.28pg.118*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-28*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-28",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGATCCAGCTGCAACAGTCTGGACCTGAGCTGGTGAAGCCTGGTGCTTCAGTGAAGATATCCTTCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGAACTGGATGAAGTAGAGCCACGGAAATAGCCTTGAGTGGATTGGATATATTGATCCTTACAATGGTGGTACTAGCTACAACCAGAAATTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCTTCCAGCACAGCCTACATGCATCTCAACAACCTGACATCTGAGGACTCTGCAGTCTATTACTTTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1261354,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1258766,
                "length": 7,
                "sequence": {
                    "nucleotide": "TACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1258759
            },
            "3 nonamer": {
                "end": 1258798,
                "length": 9,
                "sequence": {
                    "nucleotide": "TAGTAATCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1258789
            },
            "end": 1258759,
            "length": 296,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCACCTGTAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTCCTGGTTACTCATTCATTGGCTACTACATGCACTGAGTAAAGTCCTGTAAAAAGCCTTTAAAGGATTGGATATATTAATCCAGTGGTGGTACTGGCTACAATGAAAAGTTCAAGGGCAAGGCCACATTGACTGCAGACAAATCCTTCAGCACAGCCTATATGCAATTCAGCAGCCTGACATCTGAGGACTCTGTGGACTCTTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1258463
        },
        "accession": "BN000872",
        "allele": "J558.29pg.119*01",
        "confirmed": True,
        "end": 1258759,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.29pg.119",
        "id": "J558.29pg.119*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-29*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-29",
        "length": 296,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCACCTGTAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTCCTGGTTACTCATTCATTGGCTACTACATGCACTGAGTAAAGTCCTGTAAAAAGCCTTTAAAGGATTGGATATATTAATCCAGTGGTGGTACTGGCTACAATGAAAAGTTCAAGGGCAAGGCCACATTGACTGCAGACAAATCCTTCAGCACAGCCTATATGCAATTCAGCAGCCTGACATCTGAGGACTCTGTGGACTCTTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1258463,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1557786,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1557779
            },
            "3 nonamer": {
                "end": 1557818,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACT",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1557809
            },
            "end": 1557779,
            "length": 305,
            "sequence": {
                "nucleotide": "GGGTCTACTCAGAGGTTCAGCTCCAGCAGTCTGGGACTGTGCTGGCAAGGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGACTTCTGGCTACACATTTACCAGCTACTGGATGCACTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATAGGGGCTATTTATCCTGGAAATAGTGATACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCAAACTGACTGCAGTCACATCCGCCAGCACTGCCTACATGGAGCTCAGCAGCCTGACAAATGAGGACTCTGCGGTCTATTACTGTACAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1557474
        },
        "accession": "BN000872",
        "allele": "J558.3.90*01",
        "confirmed": True,
        "end": 1557779,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.3.90",
        "id": "J558.3.90*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-5*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-5",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGGTCTACTCAGAGGTTCAGCTCCAGCAGTCTGGGACTGTGCTGGCAAGGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGACTTCTGGCTACACATTTACCAGCTACTGGATGCACTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATAGGGGCTATTTATCCTGGAAATAGTGATACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCAAACTGACTGCAGTCACATCCGCCAGCACTGCCTACATGGAGCTCAGCAGCCTGACAAATGAGGACTCTGCGGTCTATTACTGTACAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1557474,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1253988,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1253981
            },
            "3 nonamer": {
                "end": 1254017,
                "length": 9,
                "sequence": {
                    "nucleotide": "TGTCAGAAA",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1254008
            },
            "end": 1253981,
            "length": 302,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGTTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGACTACTACATGCACTGGGTGAAGCAAAGTCCTGAAAAGAGTCTTGAGTGGATTGGAGAGATCAATCCTAGCACTGTTGTTACTAACTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAAGAGCCTGACATCTGAAGACTGCAGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1253678
        },
        "accession": "BN000872",
        "allele": "J558.30pg.120*01",
        "confirmed": True,
        "end": 1253981,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.30pg.120",
        "id": "J558.30pg.120*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-30*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-30",
        "length": 303,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGTTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGACTACTACATGCACTGGGTGAAGCAAAGTCCTGAAAAGAGTCTTGAGTGGATTGGAGAGATCAATCCTAGCACTGTTGTTACTAACTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAAGAGCCTGACATCTGAAGACTGCAGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1253678,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1241823,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1241816
            },
            "3 nonamer": {
                "end": 1241855,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1241846
            },
            "end": 1241816,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAATATCCTCGATTGGATTGGATATATTTATCCTTACAATGGTGTTTCTAGCTACAACCAGAAATTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCTAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1241511
        },
        "accession": "BN000872",
        "allele": "J558.31.121*01",
        "confirmed": True,
        "end": 1241816,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.31.121",
        "id": "J558.31.121*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-31*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-31",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAATATCCTCGATTGGATTGGATATATTTATCCTTACAATGGTGTTTCTAGCTACAACCAGAAATTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCTAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1241511,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1235864,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1235857
            },
            "3 nonamer": {
                "end": 1235896,
                "length": 10,
                "sequence": {
                    "nucleotide": "TCAGAAATCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1235886
            },
            "end": 1235857,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCAGCAAGACTTCTGGTTAAACTTTTACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTTGATATATTGATCCTTACGATGGTGTTACTAGCTACAATAAGAAGTTCAAGAGAAAGGCCACGTTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCTGAAACCTGACATCTGAGGAAACTGAAGTCTATTACTGTGTAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1235552
        },
        "accession": "BN000872",
        "allele": "J558.32pg.122*01",
        "confirmed": True,
        "end": 1235857,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.32pg.122",
        "id": "J558.32pg.122*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-32*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-32",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCAGCAAGACTTCTGGTTAAACTTTTACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTTGATATATTGATCCTTACGATGGTGTTACTAGCTACAATAAGAAGTTCAAGAGAAAGGCCACGTTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCTGAAACCTGACATCTGAGGAAACTGAAGTCTATTACTGTGTAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1235552,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1227436,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1227429
            },
            "3 nonamer": {
                "end": 1227467,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAAAC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1227458
            },
            "end": 1227429,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAACAATCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGTAAGGCTTCTGGATACACATTCACTGATTACTACATGCACTGGGTGAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGATATTAATCCTTACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAACAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1227125
        },
        "accession": "BN000872",
        "allele": "J558.33pg.123*01",
        "confirmed": True,
        "end": 1227429,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.33pg.123",
        "id": "J558.33pg.123*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-33*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-33",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAACAATCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGTAAGGCTTCTGGATACACATTCACTGATTACTACATGCACTGGGTGAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGAGATATTAATCCTTACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAACAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1227125,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1219897,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1219890
            },
            "3 nonamer": {
                "end": 1219928,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAATC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1219919
            },
            "end": 1219890,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAGTCTGGACCTGAGTTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACATTCACTGACTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGATATATTTATCCTAACAATGGTGGTAATGGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1219585
        },
        "accession": "BN000872",
        "allele": "J558.34.124*01",
        "confirmed": True,
        "end": 1219890,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.34.124",
        "id": "J558.34.124*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-34*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-34",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAGTCTGGACCTGAGTTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACATTCACTGACTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGATATATTTATCCTAACAATGGTGGTAATGGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAGTCCTCCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1219585,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1195781,
                "length": 7,
                "sequence": {
                    "nucleotide": "CTCAGTT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1195774
            },
            "3 nonamer": {
                "end": 1195813,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAGAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1195804
            },
            "end": 1195774,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGAAGCAGCCTGGAACTGTGGTGGTGAAACCTGGGGCTTCAGTGAAGATATCCTGCCAGGCTTCTGGTTACTCATTTACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGAAAAGAGCCTTTAGTGGATTCTACTTATTATTCCTTACAATGGTGATACTAGCAACAACCAGAAGTTCAAGGGCAAAGCAACATTGACTGTAGACAAGTCCTCCAGCACAGCCAACATGGAGCTCTGCAGCCTGACATCTGAGGACTCTACGGTCTATTACCGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1195469
        },
        "accession": "BN000872",
        "allele": "J558.35pg.125*01",
        "confirmed": True,
        "end": 1195774,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.35pg.125",
        "id": "J558.35pg.125*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-35*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-35",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGAAGCAGCCTGGAACTGTGGTGGTGAAACCTGGGGCTTCAGTGAAGATATCCTGCCAGGCTTCTGGTTACTCATTTACTGGCTACTACATGCACTGGGTGAAGCAGAGCCATGAAAAGAGCCTTTAGTGGATTCTACTTATTATTCCTTACAATGGTGATACTAGCAACAACCAGAAGTTCAAGGGCAAAGCAACATTGACTGTAGACAAGTCCTCCAGCACAGCCAACATGGAGCTCTGCAGCCTGACATCTGAGGACTCTACGGTCTATTACCGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1195469,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1191196,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1191189
            },
            "3 nonamer": {
                "end": 1191227,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAAAC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1191218
            },
            "end": 1191189,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAGTCTGGACCTGTGCTGGTGAAGCCTGGGCCTTCAGTGAAGATATCCTGTAAGGCTTCTGGATTCACATTCACTGACTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACTTGTTTATCCTTACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGGAGCTAAACAGCCTGACTTCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1190884
        },
        "accession": "BN000872",
        "allele": "J558.36.126*01",
        "confirmed": True,
        "end": 1191189,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.36.126",
        "id": "J558.36.126*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-36*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-36",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTCTCTGAGGTCCAGCTGCAACAGTCTGGACCTGTGCTGGTGAAGCCTGGGCCTTCAGTGAAGATATCCTGTAAGGCTTCTGGATTCACATTCACTGACTACTACATGCACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACTTGTTTATCCTTACAATGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGGAGCTAAACAGCCTGACTTCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1190884,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1174846,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACTGTA",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1174839
            },
            "3 nonamer": {
                "end": 1174877,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAAAC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1174868
            },
            "end": 1174839,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTGTTCTCTGAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTTACTGGCTACTTTATGAACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACGTATTAATCCTTACAATGGTGATACTTTCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCTAGCACAGCCCACATGGAGCTCCTGAGCCTGACATCTGAGGACTTTGCAGTCTATTATTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1174534
        },
        "accession": "BN000872",
        "allele": "J558.37.127*01",
        "confirmed": True,
        "end": 1174839,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.37.127",
        "id": "J558.37.127*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-37*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-37",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTGTTCTCTGAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTTACTGGCTACTTTATGAACTGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTGGATTGGACGTATTAATCCTTACAATGGTGATACTTTCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCTAGCACAGCCCACATGGAGCTCCTGAGCCTGACATCTGAGGACTTTGCAGTCTATTATTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1174534,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1161073,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1161066
            },
            "3 nonamer": {
                "end": 1161104,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAAAACT",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1161095
            },
            "end": 1161066,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTCGTGAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTAGCTACTAAATGCACTGGGTGAAGCAAAGCCATGGAAAGAGCCTTGAGTGGATTGGACTTATTATTCCTTACAATGGTGATACTGGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACTAGTCCTTCAGCACAACCTACATGGAGCTCCGCAGCCTGAAATCTGAGGACTCTGCGGTCTATTACTCTACAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1160762
        },
        "accession": "BN000872",
        "allele": "J558.38pg.128*01",
        "confirmed": True,
        "end": 1161066,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.38pg.128",
        "id": "J558.38pg.128*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-38*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-38",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTCGTGAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTAGCTACTAAATGCACTGGGTGAAGCAAAGCCATGGAAAGAGCCTTGAGTGGATTGGACTTATTATTCCTTACAATGGTGATACTGGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACTAGTCCTTCAGCACAACCTACATGGAGCTCCGCAGCCTGAAATCTGAGGACTCTGCGGTCTATTACTCTACAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1160762,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1156485,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1156478
            },
            "3 nonamer": {
                "end": 1156517,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1156508
            },
            "end": 1156478,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGTTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGCGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGACTACAACATGAACTGGGTGAAGCAGAGCAATGGAAAGAGCCTTGAGTGGATTGGAGTAATTAATCCTAACTATGGTACTACTAGCTACAATCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACCAATCTTCCAGCACAGCCTACATGCAGCTCAACAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1156173
        },
        "accession": "BN000872",
        "allele": "J558.39.129*01",
        "confirmed": True,
        "end": 1156478,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.39.129",
        "id": "J558.39.129*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-39*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-39",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGTTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGCGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGACTACAACATGAACTGGGTGAAGCAGAGCAATGGAAAGAGCCTTGAGTGGATTGGAGTAATTAATCCTAACTATGGTACTACTAGCTACAATCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACCAATCTTCCAGCACAGCCTACATGCAGCTCAACAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1156173,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1532620,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1532613
            },
            "3 nonamer": {
                "end": 1532652,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1532643
            },
            "end": 1532613,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGGGCTGAACTGGCAAAACCTGGGGCCTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTTACTAGCTACTGGATGCACTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATTGGATACATTAATCCTAGCAGTGGTTATACTAAGTACAATCAGAAGTTCAAGGACAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTGAGCAGCCTGACATATGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1532308
        },
        "accession": "BN000872",
        "allele": "J558.4.93*01",
        "confirmed": True,
        "end": 1532613,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.4.93",
        "id": "J558.4.93*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-7*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-7",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGGGCTGAACTGGCAAAACCTGGGGCCTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTTACTAGCTACTGGATGCACTGGGTAAAACAGAGGCCTGGACAGGGTCTGGAATGGATTGGATACATTAATCCTAGCAGTGGTTATACTAAGTACAATCAGAAGTTCAAGGACAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTGAGCAGCCTGACATATGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1532308,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1146143,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1146136
            },
            "3 nonamer": {
                "end": 1146175,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1146166
            },
            "end": 1146136,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAACTGCAACAGTCTGGATCTAAGGTAGTGAATGCTGGGGCTTCCATGAAGCTGTCCTGCAAATATTCTGGTTACTCATTTAGTAGATACAAAATGGAATGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTCGATTGGAGATATTAATCTTTCCAATGGTGGTACTAACTACAATGGAAAGTTCAAAAGCAAGGCCACATTGACTGTAGATATATCCTCTAGCACAGCCTATATGGAGCTAGCATATTGACATCTGAGGTCTCTGCAGTCTCTCACCATGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1145832
        },
        "accession": "BN000872",
        "allele": "J558.40pg.130*01",
        "confirmed": True,
        "end": 1146136,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.40pg.130",
        "id": "J558.40pg.130*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-40*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-40",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAACTGCAACAGTCTGGATCTAAGGTAGTGAATGCTGGGGCTTCCATGAAGCTGTCCTGCAAATATTCTGGTTACTCATTTAGTAGATACAAAATGGAATGGGTGAAGCAGAGCCATGGAAAGAGCCTTGAGTCGATTGGAGATATTAATCTTTCCAATGGTGGTACTAACTACAATGGAAAGTTCAAAAGCAAGGCCACATTGACTGTAGATATATCCTCTAGCACAGCCTATATGGAGCTAGCATATTGACATCTGAGGTCTCTGCAGTCTCTCACCATGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1145832,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1138713,
                "length": 7,
                "sequence": {
                    "nucleotide": "TACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1138706
            },
            "3 nonamer": {
                "end": 1138745,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGTAATAC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1138736
            },
            "end": 1138706,
            "length": 283,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCACCTGAAGCAGTCTGGACCTAAGGTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGTCCTGTAAAAAGCCTTCAAAGGATTGGATATATTATTCCTTACAGTGGTGGTACTGGGTACAATGAAAAGTTCAAGGGCAAGGCCACATTGTCTGCAGACAAATTCTTCAGCACAACCTATATGTACTTCAGCAGCCTGACATCTGAGGACTCTGTGGACTCTTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1138423
        },
        "accession": "BN000872",
        "allele": "J558.41pg.131*01",
        "confirmed": True,
        "end": 1138706,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.41pg.131",
        "id": "J558.41pg.131*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-41*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-41",
        "length": 283,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCACCTGAAGCAGTCTGGACCTAAGGTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGTCCTGTAAAAAGCCTTCAAAGGATTGGATATATTATTCCTTACAGTGGTGGTACTGGGTACAATGAAAAGTTCAAGGGCAAGGCCACATTGTCTGCAGACAAATTCTTCAGCACAACCTATATGTACTTCAGCAGCCTGACATCTGAGGACTCTGTGGACTCTTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1138423,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1133971,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1133964
            },
            "3 nonamer": {
                "end": 1134003,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1133994
            },
            "end": 1133964,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGAACTGGGTGAAGCAAAGTCCTGAAAAGAGCCTTGAGTGGATTGGAGAGATTAATCCTAGCACTGGTGGTACTACCTACAACCAGAAGTTCAAGGCCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAAGAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1133659
        },
        "accession": "BN000872",
        "allele": "J558.42.132*01",
        "confirmed": True,
        "end": 1133964,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.42.132",
        "id": "J558.42.132*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-42*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-42",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGAACTGGGTGAAGCAAAGTCCTGAAAAGAGCCTTGAGTGGATTGGAGAGATTAATCCTAGCACTGGTGGTACTACCTACAACCAGAAGTTCAAGGCCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAAGAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1133659,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1125134,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACGGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1125127
            },
            "3 nonamer": {
                "end": 1125166,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACA",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1125157
            },
            "end": 1125127,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCAAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAAAGTTCTGAAAAGAGCCTTGAGTGGATTGGAGAGATTAATCCTAGCACTGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTAACTGTAGACAAGTCATCCAGCACAGCCTACATGCAGCTCAAGAGCCTGACATCTGAGGACTCTGCTGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1124822
        },
        "accession": "BN000872",
        "allele": "J558.43.133*01",
        "confirmed": True,
        "end": 1125127,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.43.133",
        "id": "J558.43.133*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-43*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-43",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCAAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGTTACTCATTCACTGGCTACTACATGCACTGGGTGAAGCAAAGTTCTGAAAAGAGCCTTGAGTGGATTGGAGAGATTAATCCTAGCACTGGTGGTACTAGCTACAACCAGAAGTTCAAGGGCAAGGCCACATTAACTGTAGACAAGTCATCCAGCACAGCCTACATGCAGCTCAAGAGCCTGACATCTGAGGACTCTGCTGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1124822,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1107704,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1107697
            },
            "3 nonamer": {
                "end": 1107736,
                "length": 9,
                "sequence": {
                    "nucleotide": "AGAAACCCT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1107727
            },
            "end": 1107697,
            "length": 306,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGGGCTGAAATGGTGAGGCCTGGGGCCTCAGTGAAGATGTCCTGCAAGGCTTTTGGCTACACATTCACTGACTATGGCATGCATTGGGTGAAGCATAGTCACGGAAGGCATCCTTGAGAGGATTGGAAATATTAATACTTACTATGACAATACTGGCTACAATGAAAAGTTCAAGGGCAAGTCCAAATTGACTGTAGACAAATCCTCCAGCACAGCCTATGTGGAGTTTAGCAGAATGACATCTGAGGATTCTGTAGTCTATTACTGTGAAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1107391
        },
        "accession": "BN000872",
        "allele": "J558.44pg.134*01",
        "confirmed": True,
        "end": 1107697,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.44pg.134",
        "id": "J558.44pg.134*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-44*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-44",
        "length": 306,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGGGCTGAAATGGTGAGGCCTGGGGCCTCAGTGAAGATGTCCTGCAAGGCTTTTGGCTACACATTCACTGACTATGGCATGCATTGGGTGAAGCATAGTCACGGAAGGCATCCTTGAGAGGATTGGAAATATTAATACTTACTATGACAATACTGGCTACAATGAAAAGTTCAAGGGCAAGTCCAAATTGACTGTAGACAAATCCTCCAGCACAGCCTATGTGGAGTTTAGCAGAATGACATCTGAGGATTCTGTAGTCTATTACTGTGAAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1107391,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1101029,
                "length": 7,
                "sequence": {
                    "nucleotide": "TACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1101022
            },
            "3 nonamer": {
                "end": 1101061,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAATCT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1101052
            },
            "end": 1101022,
            "length": 306,
            "sequence": {
                "nucleotide": "TTGTTCACTCCCAGGTTCAGATGCAGCAGTCTGGGGGATTAGGGGATGAAGCCTGTGTTCTCAGTGAAGATGTCCTGTAAGGGTTCTGGCTACAGCTTTACCAACTACTATATGCCCTGGGTAAAATAGTGGACTGGACACAACCTTGAGTGGATTGGATGGATTCATCCTGGAAATGGTGATACTTACTACAATCAAAAGTTCAAGGGAAAGGCAACACTGACCAAGTACAAATTCTCCAGCACAGCCTACTTACATCACAACAGCCTGACATCTGAGCACCCAGTAGTTTATAAATATGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1100716
        },
        "accession": "BN000872",
        "allele": "J558.45pg.135*01",
        "confirmed": True,
        "end": 1101022,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.45pg.135",
        "id": "J558.45pg.135*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-45*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-45",
        "length": 306,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TTGTTCACTCCCAGGTTCAGATGCAGCAGTCTGGGGGATTAGGGGATGAAGCCTGTGTTCTCAGTGAAGATGTCCTGTAAGGGTTCTGGCTACAGCTTTACCAACTACTATATGCCCTGGGTAAAATAGTGGACTGGACACAACCTTGAGTGGATTGGATGGATTCATCCTGGAAATGGTGATACTTACTACAATCAAAAGTTCAAGGGAAAGGCAACACTGACCAAGTACAAATTCTCCAGCACAGCCTACTTACATCACAACAGCCTGACATCTGAGCACCCAGTAGTTTATAAATATGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1100716,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1093391,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1093384
            },
            "3 nonamer": {
                "end": 1093421,
                "length": 7,
                "sequence": {
                    "nucleotide": "CAAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1093414
            },
            "end": 1093384,
            "length": 305,
            "sequence": {
                "nucleotide": "GTATCCACTGCCAGGTCCAGGTGCAGCTGTCTGCAGCTGAGCTGGTGAAGCCTGGGAGTCCAGTGAAGCTGTCCTGCAAAGCTTCTGGCTACACCGTCAATGACAACTATATGGAGCAGGTAAAGCAGAGGCCTGGACAGAGCATGGAATGGATTGGATAGATTCATTTTGTATATGGTGGTACTTAATACAATGAAAAGTTCTAGGGCAAGTCCACATTAACTGTAGAAAAATCCTCCAACACAGCCTACATGGAACTCAACAGCTCGACATCTGAGGACTCTGTAGTTTATTACTGTGCATGG",
                "read frame": 1,
                "strand": True
            },
            "start": 1093079
        },
        "accession": "BN000872",
        "allele": "J558.46pg.136*01",
        "confirmed": True,
        "end": 1093384,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.46pg.136",
        "id": "J558.46pg.136*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-46*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-46",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCACTGCCAGGTCCAGGTGCAGCTGTCTGCAGCTGAGCTGGTGAAGCCTGGGAGTCCAGTGAAGCTGTCCTGCAAAGCTTCTGGCTACACCGTCAATGACAACTATATGGAGCAGGTAAAGCAGAGGCCTGGACAGAGCATGGAATGGATTGGATAGATTCATTTTGTATATGGTGGTACTTAATACAATGAAAAGTTCTAGGGCAAGTCCACATTAACTGTAGAAAAATCCTCCAACACAGCCTACATGGAACTCAACAGCTCGACATCTGAGGACTCTGTAGTTTATTACTGTGCATGG",
            "read frame": 1,
            "strand": True
        },
        "start": 1093079,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1079976,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1079969
            },
            "3 nonamer": {
                "end": 1080008,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAAAAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1079999
            },
            "end": 1079969,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTTCAGCTGCAGCAGTCTGGGGCTGAGCTAGTGAAGCCTGGAGCCTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTCACTACCTATCCTATAGAGTGGATGAAGCAGAATCATGGAAAGAGCCTAGAGTGGATTGGAAATTTTCATCCTTACAATGATGATACTAAGTACAATGAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGAAAAATCCTCTAGCACAGTCTACTTGGAGCTCAGCCGATTAACATCTGATGACTCTGCTGTTTATTACTGTGCAAGG",
                "read frame": 0,
                "strand": True
            },
            "start": 1079664
        },
        "accession": "BN000872",
        "allele": "J558.47.137*01",
        "confirmed": True,
        "end": 1079969,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.47.137",
        "id": "J558.47.137*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-47*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-47",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTTCAGCTGCAGCAGTCTGGGGCTGAGCTAGTGAAGCCTGGAGCCTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTCACTACCTATCCTATAGAGTGGATGAAGCAGAATCATGGAAAGAGCCTAGAGTGGATTGGAAATTTTCATCCTTACAATGATGATACTAAGTACAATGAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGAAAAATCCTCTAGCACAGTCTACTTGGAGCTCAGCCGATTAACATCTGATGACTCTGCTGTTTATTACTGTGCAAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1079664,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1032857,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1032850
            },
            "3 nonamer": {
                "end": 1032889,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACTC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1032880
            },
            "end": 1032850,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCTACTTCTAGGTCCAATGGCAGGAGTCAGGGACTGAGCTGGTGAGATCTGGGGCCTCAGTAATGATGTCCTGCAAGGCTTCTGGATACACATTCAGTAACTACACTATGCACTGGGTAAAGCAGAGTCATGGAAAAGGCCATGAAAGGATTGGATATATTGATAATTACTATTGTAGCACTGACTACAGTGAAAAGTTCAAGATCAAGGCCACATTGACTGTAAACAAATCCTGCAGAACAGCCTATGTCAAGCTCAGCAGACTGACATCTGAGGACTCTGCAGTCTATTATTGTGTAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1032545
        },
        "accession": "BN000872",
        "allele": "J558.48pg.140*01",
        "confirmed": True,
        "end": 1032850,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.48pg.140",
        "id": "J558.48pg.140*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-48*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-48",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCTACTTCTAGGTCCAATGGCAGGAGTCAGGGACTGAGCTGGTGAGATCTGGGGCCTCAGTAATGATGTCCTGCAAGGCTTCTGGATACACATTCAGTAACTACACTATGCACTGGGTAAAGCAGAGTCATGGAAAAGGCCATGAAAGGATTGGATATATTGATAATTACTATTGTAGCACTGACTACAGTGAAAAGTTCAAGATCAAGGCCACATTGACTGTAAACAAATCCTGCAGAACAGCCTATGTCAAGCTCAGCAGACTGACATCTGAGGACTCTGCAGTCTATTATTGTGTAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1032545,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1015861,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1015854
            },
            "3 nonamer": {
                "end": 1015893,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACTC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1015884
            },
            "end": 1015854,
            "length": 305,
            "sequence": {
                "nucleotide": "GAGTCCACTCACAGCGTGAGCTGCAGCAGTCTGGAGCTGAGTTGGTGAGACCTGGGTCCTCAGTGAAGTTGTCCTGCAAGGATTCTTACTTTGCCTTCATGGCCAGTGCTATGCACTGGGTGAAGCAAAGACCTGGACATGGCCTGGAATGGATAGGATCTTTTACTATGTACAGTGATGCTACTGAGTACAGTGAAAACTTCAAGGGCAAGGCCACATTGACTGCAAACACATCCTCGAGCACAGCCTACATGGAACTCAGCAGCCTCACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1015549
        },
        "accession": "BN000872",
        "allele": "J558.49.141*01",
        "confirmed": True,
        "end": 1015854,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.49.141",
        "id": "J558.49.141*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-49*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-49",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGTCCACTCACAGCGTGAGCTGCAGCAGTCTGGAGCTGAGTTGGTGAGACCTGGGTCCTCAGTGAAGTTGTCCTGCAAGGATTCTTACTTTGCCTTCATGGCCAGTGCTATGCACTGGGTGAAGCAAAGACCTGGACATGGCCTGGAATGGATAGGATCTTTTACTATGTACAGTGATGCTACTGAGTACAGTGAAAACTTCAAGGGCAAGGCCACATTGACTGCAAACACATCCTCGAGCACAGCCTACATGGAACTCAGCAGCCTCACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1015549,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 951337,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 951330
            },
            "3 nonamer": {
                "end": 951369,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 951360
            },
            "end": 951330,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCAGTGGGTAAAACAGAGGCCTGGACAGGGCCTTGAGTGGATCGGAGAGATTGATCCTTCTGATAGCTATACTAACTACAATCAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 951025
        },
        "accession": "BN000872",
        "allele": "J558.50.143*01",
        "confirmed": True,
        "end": 951330,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.50.143",
        "id": "J558.50.143*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-50*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-50",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCAGTGGGTAAAACAGAGGCCTGGACAGGGCCTTGAGTGGATCGGAGAGATTGATCCTTCTGATAGCTATACTAACTACAATCAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 951025,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 939716,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 939709
            },
            "3 nonamer": {
                "end": 939748,
                "length": 9,
                "sequence": {
                    "nucleotide": "CATAAACTC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 939739
            },
            "end": 939709,
            "length": 311,
            "sequence": {
                "nucleotide": "GTATCTACTCCCAGGTCCAGCTTCCTCAGTCTGGTTCTGAGGTGGGGCGGACTGGTGCCTCAGTGAAGATGTTCTGCAAGGCTCCTGGCTACACATTCACTAACTACTATATGTATTGATTAAAGCAGAGTCATGGAGATAGCCTAGAGTGGATTTGATATATTTATCCTGGAAATGGTCTTACTAGCTATGCCAAGAAGTTCAAAGGCAAGGCCACATTGACTATAGACAATTCAGCCAGCACAGCCTACATGCAGCTCAGCAGCATGACATCTGAAGCCTCTGATGACTATTGTTGTGCTAGACAAGTG",
                "read frame": 1,
                "strand": True
            },
            "start": 939398
        },
        "accession": "BN000872",
        "allele": "J558.51pg.144*01",
        "confirmed": True,
        "end": 939709,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.51pg.144",
        "id": "J558.51pg.144*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-51*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-51",
        "length": 311,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCTACTCCCAGGTCCAGCTTCCTCAGTCTGGTTCTGAGGTGGGGCGGACTGGTGCCTCAGTGAAGATGTTCTGCAAGGCTCCTGGCTACACATTCACTAACTACTATATGTATTGATTAAAGCAGAGTCATGGAGATAGCCTAGAGTGGATTTGATATATTTATCCTGGAAATGGTCTTACTAGCTATGCCAAGAAGTTCAAAGGCAAGGCCACATTGACTATAGACAATTCAGCCAGCACAGCCTACATGCAGCTCAGCAGCATGACATCTGAAGCCTCTGATGACTATTGTTGTGCTAGACAAGTG",
            "read frame": 1,
            "strand": True
        },
        "start": 939398,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 925598,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 925591
            },
            "3 nonamer": {
                "end": 925630,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 925621
            },
            "end": 925591,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGTCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCATTGGGTGAAGCAGAGGCCTATACAAGGCCTTGAATGGATTGGTAACATTGACCCTTCTGATAGTGAAACTCACTACAATCAAAAGTTCAAGGACAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 925286
        },
        "accession": "BN000872",
        "allele": "J558.52.145*01",
        "confirmed": True,
        "end": 925591,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.52.145",
        "id": "J558.52.145*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-52*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-52",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGTCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCATTGGGTGAAGCAGAGGCCTATACAAGGCCTTGAATGGATTGGTAACATTGACCCTTCTGATAGTGAAACTCACTACAATCAAAAGTTCAAGGACAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 925286,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 912682,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 912675
            },
            "3 nonamer": {
                "end": 912714,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 912705
            },
            "end": 912675,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAAATATTAATCCTAGCAATGGTGGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 912370
        },
        "accession": "BN000872",
        "allele": "J558.53.146*01",
        "confirmed": True,
        "end": 912675,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.53.146",
        "id": "J558.53.146*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-53*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-53",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAAATATTAATCCTAGCAATGGTGGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 912370,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 877410,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 877403
            },
            "3 nonamer": {
                "end": 877442,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 877433
            },
            "end": 877403,
            "length": 294,
            "sequence": {
                "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTAAGGCCTGGGACTTCAGTGAAGGTGTCCTGCAAGGCTTCTGGATACGCCTTCACTAATTACTTGATAGAGTGGGTAAAGCAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAGTGATTAATCCTGGAAGTGGTGGTACTAACTACAATGAGAAGTTCAAGGGCAAGGCAACACTGACTGCAGACAAATCCTCCAGCACTGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 877109
        },
        "accession": "BN000872",
        "allele": "J558.54.148*01",
        "confirmed": True,
        "end": 877403,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.54.148",
        "id": "J558.54.148*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-54*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-54",
        "length": 294,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTAAGGCCTGGGACTTCAGTGAAGGTGTCCTGCAAGGCTTCTGGATACGCCTTCACTAATTACTTGATAGAGTGGGTAAAGCAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAGTGATTAATCCTGGAAGTGGTGGTACTAACTACAATGAGAAGTTCAAGGGCAAGGCAACACTGACTGCAGACAAATCCTCCAGCACTGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 877109,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 862990,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 862983
            },
            "3 nonamer": {
                "end": 863022,
                "length": 9,
                "sequence": {
                    "nucleotide": "CCAAAACCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 863013
            },
            "end": 862983,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATAACCTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAGATATTTATCCTGGTAGTGGTAGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 862678
        },
        "accession": "BN000872",
        "allele": "J558.55.149*01",
        "confirmed": True,
        "end": 862983,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.55.149",
        "id": "J558.55.149*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-55*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-55",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATAACCTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAGATATTTATCCTGGTAGTGGTAGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 862678,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 828283,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 828276
            },
            "3 nonamer": {
                "end": 828315,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACAC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 828306
            },
            "end": 828276,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTGCCAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTCCTGGCTATACCTTCACCAGCCACTGGATGCAGTGGGTAAGACAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAGAGATTTTTCCTGGAAGTGGTAGTACTTATTATAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 827971
        },
        "accession": "BN000872",
        "allele": "J558.56.150*01",
        "confirmed": True,
        "end": 828276,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.56.150",
        "id": "J558.56.150*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-56*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-56",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTGCCAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTCCTGGCTATACCTTCACCAGCCACTGGATGCAGTGGGTAAGACAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAGAGATTTTTCCTGGAAGTGGTAGTACTTATTATAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 827971,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 802661,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGGT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 802654
            },
            "3 nonamer": {
                "end": 802693,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACA",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 802684
            },
            "end": 802654,
            "length": 303,
            "sequence": {
                "nucleotide": "ATATCCACTCCCAGGTCTAGCTGCATCAATCTGGGGCTGAGCTGGTGAAGCCTGGGATCTCTGTGAAGATGTTATACATGGTTTCTGGTTATACATTTACTGAATACTACATGCCCTGGTCAAGCAGAATCATGGAAAGATCCTTGAGTGGATTGGAAATGTTAATATTTACAATTGTGGTATAACTACAATAAAAATTTCAAGGACAAGGACACATCAACTGTAGACTATTCCTCCAGTACAGCCTATATGTTGCTTGGCAAAGTGACATCTGAGGATTCTAAGGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 802351
        },
        "accession": "BN000872",
        "allele": "J558.57pg.152*01",
        "confirmed": True,
        "end": 802654,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.57pg.152",
        "id": "J558.57pg.152*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-57*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-57",
        "length": 303,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATATCCACTCCCAGGTCTAGCTGCATCAATCTGGGGCTGAGCTGGTGAAGCCTGGGATCTCTGTGAAGATGTTATACATGGTTTCTGGTTATACATTTACTGAATACTACATGCCCTGGTCAAGCAGAATCATGGAAAGATCCTTGAGTGGATTGGAAATGTTAATATTTACAATTGTGGTATAACTACAATAAAAATTTCAAGGACAAGGACACATCAACTGTAGACTATTCCTCCAGTACAGCCTATATGTTGCTTGGCAAAGTGACATCTGAGGATTCTAAGGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 802351,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 758919,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 758912
            },
            "3 nonamer": {
                "end": 758951,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACAC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 758942
            },
            "end": 758912,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTGAGGTCCAGCTTCAGCAGTCTGGAGCTGAGCTGGTGAGGCCTGGGTCCTCAGTGAAGATGTCCTGCAAGACTTCTGGATATACATTCACAAGCTACGGTATAAACTGGGTGAAGCAGAGGCCTGGACAGGGCCTGGAATGGATTGGATATATTTATATTGGAAATGGTTATACTGAGTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTTCAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCAATCTATTTCTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 758607
        },
        "accession": "BN000872",
        "allele": "J558.58.154*01",
        "confirmed": True,
        "end": 758912,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.58.154",
        "id": "J558.58.154*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-58*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-58",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTGAGGTCCAGCTTCAGCAGTCTGGAGCTGAGCTGGTGAGGCCTGGGTCCTCAGTGAAGATGTCCTGCAAGACTTCTGGATATACATTCACAAGCTACGGTATAAACTGGGTGAAGCAGAGGCCTGGACAGGGCCTGGAATGGATTGGATATATTTATATTGGAAATGGTTATACTGAGTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTTCAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCAATCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 758607,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 736003,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 735996
            },
            "3 nonamer": {
                "end": 736035,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 736026
            },
            "end": 735996,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGACTTCAGTGAAGTTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTAAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATCGGAGTGATTGATCCTTCTGATAGTTATACTAACTACAATCAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 735691
        },
        "accession": "BN000872",
        "allele": "J558.59.155*01",
        "confirmed": True,
        "end": 735996,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.59.155",
        "id": "J558.59.155*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-59*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-59",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGACTTCAGTGAAGTTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTAAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATCGGAGTGATTGATCCTTCTGATAGTTATACTAACTACAATCAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 735691,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1515506,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACTGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1515499
            },
            "3 nonamer": {
                "end": 1515538,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1515529
            },
            "end": 1515499,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCACAGGATCAGCTGCATCAGTCTGAAGCTGAGCTGCAGCAACCTGGGACATCAGTGAAGATGCCCTGAAAGGCTACTGGCTACACCTTCACTAAGTATCGAATGTGTTGGGTGAGGCAGAAGCTTGGACAGGGCCTGGAATGGATTGCATCTGTTGATCCTGGAAATAGTAATACTGAATACAATCAGAAGTTCAAAGGCAAGGCCACACTAACAGAACACAAATCCTCCAGCACAGCCTACATAGAGCTTAGCAACCTGACCTCTGAGGACTCTGCTGTCTATTACTGTACAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1515194
        },
        "accession": "BN000872",
        "allele": "J558.5pg.94*01",
        "confirmed": True,
        "end": 1515499,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.5pg.94",
        "id": "J558.5pg.94*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-8*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-8",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCACAGGATCAGCTGCATCAGTCTGAAGCTGAGCTGCAGCAACCTGGGACATCAGTGAAGATGCCCTGAAAGGCTACTGGCTACACCTTCACTAAGTATCGAATGTGTTGGGTGAGGCAGAAGCTTGGACAGGGCCTGGAATGGATTGCATCTGTTGATCCTGGAAATAGTAATACTGAATACAATCAGAAGTTCAAAGGCAAGGCCACACTAACAGAACACAAATCCTCCAGCACAGCCTACATAGAGCTTAGCAACCTGACCTCTGAGGACTCTGCTGTCTATTACTGTACAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1515194,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1487546,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1487539
            },
            "3 nonamer": {
                "end": 1487578,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACTC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1487569
            },
            "end": 1487539,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTTCAGCTGCAGCAGTCTGGAGCTGAGCTGATGAAGCCTGGGGCCTCAGTGAAGCTTTCCTGCAAGGCTACTGGCTACACATTCACTGGCTACTGGATAGAGTGGGTAAAGCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGAGATTTTACCTGGAAGTGGTAGTACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACATTCACTGCAGATACATCCTCCAACACAGCCTACATGCAACTCAGCAGCCTGACAACTGAGGACTCTGCCATCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1487234
        },
        "accession": "BN000872",
        "allele": "J558.6.96*01",
        "confirmed": True,
        "end": 1487539,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.6.96",
        "id": "J558.6.96*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-9*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-9",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTTCAGCTGCAGCAGTCTGGAGCTGAGCTGATGAAGCCTGGGGCCTCAGTGAAGCTTTCCTGCAAGGCTACTGGCTACACATTCACTGGCTACTGGATAGAGTGGGTAAAGCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGAGATTTTACCTGGAAGTGGTAGTACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACATTCACTGCAGATACATCCTCCAACACAGCCTACATGCAACTCAGCAGCCTGACAACTGAGGACTCTGCCATCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1487234,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 726554,
                "length": 7,
                "sequence": {
                    "nucleotide": "CATAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 726547
            },
            "3 nonamer": {
                "end": 726586,
                "length": 9,
                "sequence": {
                    "nucleotide": "CATAAACTC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 726577
            },
            "end": 726547,
            "length": 311,
            "sequence": {
                "nucleotide": "GTATCTACTCCCAGGTCCAGCTTCCTCAGTCTGGATCTGAGGTGGGGCGGACTGGTGCCTCAGTGAAGATGTTCTGCAAGGCTTCTGGCTACACATTCACTAACTACTATATGCATTGATTAAAGCAGAGTCATGGAGATAGCCTAGAGTGGATTTGATATATTTATCCTGGAAATGGTCTTACTAGCTATGCCAAGAAGTTCAAAGGCAAGGCCACATTGACTATAGACAATTCAGCCAGCACAGCCTACATGCAGCTCAGCAGCATGACATCTGAAGCCTCTGATGACTATTACTGTGCTAGACAAGTG",
                "read frame": 1,
                "strand": True
            },
            "start": 726236
        },
        "accession": "BN000872",
        "allele": "J558.60pg.156*01",
        "confirmed": True,
        "end": 726547,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.60pg.156",
        "id": "J558.60pg.156*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-60*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-60",
        "length": 311,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCTACTCCCAGGTCCAGCTTCCTCAGTCTGGATCTGAGGTGGGGCGGACTGGTGCCTCAGTGAAGATGTTCTGCAAGGCTTCTGGCTACACATTCACTAACTACTATATGCATTGATTAAAGCAGAGTCATGGAGATAGCCTAGAGTGGATTTGATATATTTATCCTGGAAATGGTCTTACTAGCTATGCCAAGAAGTTCAAAGGCAAGGCCACATTGACTATAGACAATTCAGCCAGCACAGCCTACATGCAGCTCAGCAGCATGACATCTGAAGCCTCTGATGACTATTACTGTGCTAGACAAGTG",
            "read frame": 1,
            "strand": True
        },
        "start": 726236,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 711945,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGCG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 711938
            },
            "3 nonamer": {
                "end": 711977,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 711968
            },
            "end": 711938,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGTCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGGATTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAATGGATTGGTAACATTTACCCTTCTGATAGTGAAACTCACTACAATCAAAAGTTCAAGGACAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 711633
        },
        "accession": "BN000872",
        "allele": "J558.61.157*01",
        "confirmed": True,
        "end": 711938,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.61.157",
        "id": "J558.61.157*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-61*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-61",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGTCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGGATTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAATGGATTGGTAACATTTACCCTTCTGATAGTGAAACTCACTACAATCAAAAGTTCAAGGACAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 711633,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 699099,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 699092
            },
            "3 nonamer": {
                "end": 699131,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 699122
            },
            "end": 699092,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCGCAAGGCTTCTGGCTACACTTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGAGAAGGCCTTGAGTGGATTGGAAATATTTATCCTGGTAGTAGTAGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 698788
        },
        "accession": "BN000872",
        "allele": "J558.62pg.158*01",
        "confirmed": True,
        "end": 699092,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.62pg.158",
        "id": "J558.62pg.158*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-62*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-62",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCGCAAGGCTTCTGGCTACACTTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGAGAAGGCCTTGAGTGGATTGGAAATATTTATCCTGGTAGTAGTAGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 698788,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 684389,
                "length": 6,
                "sequence": {
                    "nucleotide": "CAGTGT",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 684383
            },
            "3 nonamer": {
                "end": 684421,
                "length": 9,
                "sequence": {
                    "nucleotide": "AGAAACCCT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 684412
            },
            "end": 684383,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTGCCAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCAGTGGGTAAAACAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAGAGATTTTTCCTGGAAGTGGTAGTACTTATTATAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACAGCAGAGAACTCTGCAATCTATCTGTGCAAGGAA",
                "read frame": 0,
                "strand": True
            },
            "start": 684078
        },
        "accession": "BN000872",
        "allele": "J558.63pg.159*01",
        "confirmed": True,
        "end": 684383,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.63pg.159",
        "id": "J558.63pg.159*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-62-1*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-62-1",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTGCCAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCAGTGGGTAAAACAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAGAGATTTTTCCTGGAAGTGGTAGTACTTATTATAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACAGCAGAGAACTCTGCAATCTATCTGTGCAAGGAA",
            "read frame": 0,
            "strand": True
        },
        "start": 684078,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 624669,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 624662
            },
            "3 nonamer": {
                "end": 624701,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 624692
            },
            "end": 624662,
            "length": 313,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTGAAACCCGGGGCATCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACTGAGTATACTATACACTGGGTAAAGCAGAGGTCTGGACAGGGTCTTGAGTGGATTGGGTGGTTTTACCCTGGAAGTGGTAGTATAAAGTACAATGAGAAATTCAAGGACAAGGCCACATTGACTGCGGACAAATCCTCCAGCACAGTCTATATGGAGCTTAGTAGATTGACATCTGAAGACTCTGCGGTCTATTTCTGTGCAAGACACGAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 624349
        },
        "accession": "BN000872",
        "allele": "J558.64.162*01",
        "confirmed": True,
        "end": 624662,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.64.162",
        "id": "J558.64.162*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-62-2*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-62-2",
        "length": 313,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTGAAACCCGGGGCATCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACTGAGTATACTATACACTGGGTAAAGCAGAGGTCTGGACAGGGTCTTGAGTGGATTGGGTGGTTTTACCCTGGAAGTGGTAGTATAAAGTACAATGAGAAATTCAAGGACAAGGCCACATTGACTGCGGACAAATCCTCCAGCACAGTCTATATGGAGCTTAGTAGATTGACATCTGAAGACTCTGCGGTCTATTTCTGTGCAAGACACGAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 624349,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 610079,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 610072
            },
            "3 nonamer": {
                "end": 610111,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 610102
            },
            "end": 610072,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTTCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCAGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACGAGGCCTTGAGTGGATTGGAAGGATTGATCCTAATAGTGGTGGTACTAAGTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAACCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 609767
        },
        "accession": "BN000872",
        "allele": "J558.65.163*01",
        "confirmed": True,
        "end": 610072,
        "family": "J558",
        "framed": True,
        "functionality": "O",
        "gene": "J558.65.163",
        "id": "J558.65.163*01-C57BL/6",
        "identified": True,
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-62-3",
        "length": 305,
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTTCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCAGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACGAGGCCTTGAGTGGATTGGAAGGATTGATCCTAATAGTGGTGGTACTAAGTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAACCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 609767,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 575452,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 575445
            },
            "3 nonamer": {
                "end": 575484,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 575475
            },
            "end": 575445,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTAAGGCCTGGGACTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACCTTCACTAACTACTGGATAGGTTGGGCAAAGCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGATATTTACCCTGGAGGTGGTTATACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGTTCAGCAGCCTGACATCTGAGGACTCTGCCATCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 575140
        },
        "accession": "BN000872",
        "allele": "J558.66.165*01",
        "confirmed": True,
        "end": 575445,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.66.165",
        "id": "J558.66.165*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-63*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-63",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTAAGGCCTGGGACTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACCTTCACTAACTACTGGATAGGTTGGGCAAAGCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGATATTTACCCTGGAGGTGGTTATACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGTTCAGCAGCCTGACATCTGAGGACTCTGCCATCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 575140,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 563540,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 563533
            },
            "3 nonamer": {
                "end": 563572,
                "length": 9,
                "sequence": {
                    "nucleotide": "CTGAAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 563563
            },
            "end": 563533,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGTTGTCCTGCAAGGCTTCTGGCTACACTTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAATGATTCATCCTAATAGTGGTAGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAACTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 563228
        },
        "accession": "BN000872",
        "allele": "J558.67.166*01",
        "confirmed": True,
        "end": 563533,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.67.166",
        "id": "J558.67.166*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-64*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-64",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGTTGTCCTGCAAGGCTTCTGGCTACACTTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAATGATTCATCCTAATAGTGGTAGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAACTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 563228,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 539029,
                "length": 7,
                "sequence": {
                    "nucleotide": "TAAGATG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 539022
            },
            "3 nonamer": {
                "end": 539061,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAACCCT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 539052
            },
            "end": 539022,
            "length": 303,
            "sequence": {
                "nucleotide": "GTATTCACTCCCAGGTCTAGCTGCAGCAGTCAGGGGCTGAGCTGGTGAAGCCTGATGGCTCTGTGAAGATGTCATGCAAGGCTTCTGGTTATACATTTACTGAATACCACATGCTCTAGTCAAGCAGAATCATGGAAAGACCCTTGAATGGATTTTCAATATTAATACTTAAAATGGTGGTATAACTACAGTGAAAATTTCAAGGGCAAGGGTACATTAACTGTAGACAAATCCTCCAGCACACCCTATTTGTTGCTTAGCAAATTGACATCTGAGGATTCTGTGGTCTATTACTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 538719
        },
        "accession": "BN000872",
        "allele": "J558.68pg.168*01",
        "confirmed": True,
        "end": 539022,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.68pg.168",
        "id": "J558.68pg.168*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-65*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-65",
        "length": 303,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATTCACTCCCAGGTCTAGCTGCAGCAGTCAGGGGCTGAGCTGGTGAAGCCTGATGGCTCTGTGAAGATGTCATGCAAGGCTTCTGGTTATACATTTACTGAATACCACATGCTCTAGTCAAGCAGAATCATGGAAAGACCCTTGAATGGATTTTCAATATTAATACTTAAAATGGTGGTATAACTACAGTGAAAATTTCAAGGGCAAGGGTACATTAACTGTAGACAAATCCTCCAGCACACCCTATTTGTTGCTTAGCAAATTGACATCTGAGGATTCTGTGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 538719,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 477975,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 477968
            },
            "3 nonamer": {
                "end": 478007,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 477998
            },
            "end": 477968,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCATTGCCAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACAGCTTCACAAGCTACTATATACACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTATCCTGGAAGTGGTAATACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACGGCAGACACATCCTCCAGCACTGCCTACATGCAGCTCAGCAGCCTAACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 477663
        },
        "accession": "BN000872",
        "allele": "J558.69.170*01",
        "confirmed": True,
        "end": 477968,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.69.170",
        "id": "J558.69.170*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-66*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-66",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCATTGCCAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACAGCTTCACAAGCTACTATATACACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTATCCTGGAAGTGGTAATACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACGGCAGACACATCCTCCAGCACTGCCTACATGCAGCTCAGCAGCCTAACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 477663,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 467145,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 467138
            },
            "3 nonamer": {
                "end": 467177,
                "length": 9,
                "sequence": {
                    "nucleotide": "CACAAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 467168
            },
            "end": 467138,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTGCACTCCCAGGTCCAGCTGCAGCAGTCTGGGCCTGAGCTGGTGAGGCCTGGGGTCTCAGTGAAGATTTCCTGCAAGGGTTCCGGCTACACATTCACTGATTATGCTATGCACTGGGTGAAACAGAGTCATGCAAAGAGTCTAGAGTGGATTGGAGTTATTAGTACTTACTATGGTGATGCTAGCTACAACCAGAAGTTCAAGGACAAGGCCACAATGACTGTAGACAAATCCTCCAGCACAGCCTATATGGAACTTGCCAGACTGACATCTGAGGACTCTGCCGTCTATTACTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 466833
        },
        "accession": "BN000872",
        "allele": "J558.70pg.171*01",
        "confirmed": True,
        "end": 467138,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.70pg.171",
        "id": "J558.70pg.171*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-67*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-67",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTGCACTCCCAGGTCCAGCTGCAGCAGTCTGGGCCTGAGCTGGTGAGGCCTGGGGTCTCAGTGAAGATTTCCTGCAAGGGTTCCGGCTACACATTCACTGATTATGCTATGCACTGGGTGAAACAGAGTCATGCAAAGAGTCTAGAGTGGATTGGAGTTATTAGTACTTACTATGGTGATGCTAGCTACAACCAGAAGTTCAAGGACAAGGCCACAATGACTGTAGACAAATCCTCCAGCACAGCCTATATGGAACTTGCCAGACTGACATCTGAGGACTCTGCCGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 466833,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 459737,
            "length": 294,
            "sequence": {
                "nucleotide": "ATGTTAACTCCCTGGTTCAGCTGCAGCAATGTGAAGCTGAGGTTGTGAGACTCAGTGATGGTGTCCTGCTCAGCTTCTGTCTACATATTTAGCAACTACTATGAAGTGTATAAAGCAGAGGCCTGGACAGGCTCTTGAGTGGATTCGATGAATTGATTCTTGAAGTTTTGGTACTACCTATAATCAGAAGTTCAAAGGCCAGGCCACATCTACTGTAGACAAATCCTCCATCACAGCCTACCTGCAACTAATATCCTGACATCTGAGAACTCTGAAGACTATTACTGTGAAGAC",
                "read frame": 1,
                "strand": True
            },
            "start": 459443
        },
        "accession": "BN000872",
        "allele": "J558.71pg.172*01",
        "confirmed": True,
        "end": 459738,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.71pg.172",
        "id": "J558.71pg.172*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-68*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-68",
        "length": 284,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CTGGTTCAGCTGCAGCAATGTGAAGCTGAGGTTGTGAGACTCAGTGATGGTGTCCTGCTCAGCTTCTGTCTACATATTTAGCAACTACTATGAAGTGTATAAAGCAGAGGCCTGGACAGGCTCTTGAGTGGATTCGATGAATTGATTCTTGAAGTTTTGGTACTACCTATAATCAGAAGTTCAAAGGCCAGGCCACATCTACTGTAGACAAATCCTCCATCACAGCCTACCTGCAACTAATATCCTGACATCTGAGAACTCTGAAGACTATTACTGTGAAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 459454,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 447924,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 447917
            },
            "3 nonamer": {
                "end": 447956,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 447947
            },
            "end": 447917,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGATGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATCGGAGAGATTGATCCTTCTGATAGTTATACTAACTACAATCAAAAGTTCAAGGGCAAGTCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 447612
        },
        "accession": "BN000872",
        "allele": "J558.72.173*01",
        "confirmed": True,
        "end": 447917,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.72.173",
        "id": "J558.72.173*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-69*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-69",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGATGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATCGGAGAGATTGATCCTTCTGATAGTTATACTAACTACAATCAAAAGTTCAAGGGCAAGTCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 447612,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 379196,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 379189
            },
            "3 nonamer": {
                "end": 379228,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 379219
            },
            "end": 379189,
            "length": 305,
            "sequence": {
                "nucleotide": "CTGTCCACTGCCAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAGGCCTGGGACTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTATACCTTCTTCACCTACTGGATGAACTGGGTGTAGCAGAGGCCTGGACAGGGCCTTGAGTGGATTGGACAGATTTTTCCTGCAAGTGGTAGTACTTACTACAATGAGATGTACAAGGACAAGGCCGCATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAAGACACTGCTGTCTATTTCTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 378884
        },
        "accession": "BN000872",
        "allele": "J558.73pg.175*01",
        "confirmed": True,
        "end": 379189,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.73pg.175",
        "id": "J558.73pg.175*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-70*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-70",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CTGTCCACTGCCAGGTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAGGCCTGGGACTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTATACCTTCTTCACCTACTGGATGAACTGGGTGTAGCAGAGGCCTGGACAGGGCCTTGAGTGGATTGGACAGATTTTTCCTGCAAGTGGTAGTACTTACTACAATGAGATGTACAAGGACAAGGCCGCATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAAGACACTGCTGTCTATTTCTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 378884,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 328880,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 328873
            },
            "3 nonamer": {
                "end": 328912,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 328903
            },
            "end": 328873,
            "length": 313,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTGAAACCCGGGGCATCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACTGAGTATACTATACACTGGGTAAAGCAGAGGTCTGGACAGGGTCTTGAGTGGATTGGGTGGTTTTACCCTGGAAGTGGTAGTATAAAGTACAATGAGAAATTCAAGGACAAGGCCACATTGACTGCGGACAAATCCTCCAGCACAGTCTATATGGAGCTTAGTAGATTGACATCTGAAGACTCTGCGGTCTATTTCTGTGCAAGACACGAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 328560
        },
        "accession": "BN000872",
        "allele": "J558.74.176*01",
        "confirmed": True,
        "end": 328873,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.74.176",
        "id": "J558.74.176*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-71*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-71",
        "length": 313,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTGAAACCCGGGGCATCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACTGAGTATACTATACACTGGGTAAAGCAGAGGTCTGGACAGGGTCTTGAGTGGATTGGGTGGTTTTACCCTGGAAGTGGTAGTATAAAGTACAATGAGAAATTCAAGGACAAGGCCACATTGACTGCGGACAAATCCTCCAGCACAGTCTATATGGAGCTTAGTAGATTGACATCTGAAGACTCTGCGGTCTATTTCTGTGCAAGACACGAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 328560,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 313101,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 313094
            },
            "3 nonamer": {
                "end": 313133,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 313124
            },
            "end": 313094,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACGAGGCCTTGAGTGGATTGGAAGGATTGATCCTAATAGTGGTGGTACTAAGTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAACCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 312789
        },
        "accession": "BN000872",
        "allele": "J558.75.177*01",
        "confirmed": True,
        "end": 313094,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.75.177",
        "id": "J558.75.177*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-72*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-72",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACGAGGCCTTGAGTGGATTGGAAGGATTGATCCTAATAGTGGTGGTACTAAGTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAACCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 312789,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 281471,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 281464
            },
            "3 nonamer": {
                "end": 281503,
                "length": 8,
                "sequence": {
                    "nucleotide": "AGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 281495
            },
            "end": 281464,
            "length": 306,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTAAGGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTCACTAACTACTGGATAGGTTGGGTAAAGCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGATATTTACCCTGGAGATGGTTATACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGTGCATGCAGCCTGACCTCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 281158
        },
        "accession": "BN000872",
        "allele": "J558.76pg.179*01",
        "confirmed": True,
        "end": 281464,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.76pg.179",
        "id": "J558.76pg.179*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-73*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-73",
        "length": 306,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAGCTGCAGCAGTCTGGAGCTGAGCTGGTAAGGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACCTTCACTAACTACTGGATAGGTTGGGTAAAGCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGATATTTACCCTGGAGATGGTTATACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGTGCATGCAGCCTGACCTCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 281158,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 269414,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 269407
            },
            "3 nonamer": {
                "end": 269446,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 269437
            },
            "end": 269407,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGGTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGCCAAGGCCTTGAGTGGATTGGAAGGATTCATCCTTCTGATAGTGATACTAACTACAATCAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAATA",
                "read frame": 1,
                "strand": True
            },
            "start": 269102
        },
        "accession": "BN000872",
        "allele": "J558.77.180*01",
        "confirmed": True,
        "end": 269407,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.77.180",
        "id": "J558.77.180*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-74*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-74",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTCCAACTGCAGCAGCCTGGGGCTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGGTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGCCAAGGCCTTGAGTGGATTGGAAGGATTCATCCTTCTGATAGTGATACTAACTACAATCAAAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAATA",
            "read frame": 1,
            "strand": True
        },
        "start": 269102,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 238111,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 238104
            },
            "3 nonamer": {
                "end": 238143,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 238134
            },
            "end": 238104,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCATTGCCAGGTCCAGCTACAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACACCTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTTTCCTGGAAGTGGTAGTACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTTACTGTAGACAAATCCTCCAGCACAGCCTACATGTTGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 237799
        },
        "accession": "BN000872",
        "allele": "J558.78.182*01",
        "confirmed": True,
        "end": 238104,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.78.182",
        "id": "J558.78.182*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-75*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-75",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCATTGCCAGGTCCAGCTACAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACACCTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTTTCCTGGAAGTGGTAGTACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTTACTGTAGACAAATCCTCCAGCACAGCCTACATGTTGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 237799,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 224178,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 224171
            },
            "3 nonamer": {
                "end": 224210,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 224201
            },
            "end": 224171,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTGTCAGGTCCAGCTGAAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACTTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGCAAGGATTTATCCTGGAAGTGGTAATACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGAAAAATCCTCCAGCACTGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCTGTCTATTTCTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 223866
        },
        "accession": "BN000872",
        "allele": "J558.79.184*01",
        "confirmed": True,
        "end": 224171,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.79.184",
        "id": "J558.79.184*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-76*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-76",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTGTCAGGTCCAGCTGAAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACTTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGCAAGGATTTATCCTGGAAGTGGTAATACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGAAAAATCCTCCAGCACTGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCTGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 223866,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1473236,
                "length": 7,
                "sequence": {
                    "nucleotide": "TACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1473229
            },
            "3 nonamer": {
                "end": 1473268,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACTC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1473259
            },
            "end": 1473229,
            "length": 318,
            "sequence": {
                "nucleotide": "TTGTCCACTCCCCAGCCCAGCTGCAGCATTCTGGGGATGACATGGAGGAGCCTGGGTCCTCAGTGAAGTTTTCCTGCATGGCTTCTGGCTATACCTTCACTGACCACTCTATGCATGTGGTAAAACAATGAAACAGAGGCCTGGACAGGGACTGGAGTGGATTGAATGGGTTGGACCTACATGTGGTGGTACTGTATATGCTAGGAAGTTTCAAGGCAAGGCCACACTGACTGTAGACAAATCAGCCATCACAACCTACATGCAGCTCAGTAGCCTGACATCTGAAAACTCTGCAGTCTATTACTGTGCCATGCAAGGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1472910
        },
        "accession": "BN000872",
        "allele": "J558.7pg.97*01",
        "confirmed": True,
        "end": 1473229,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.7pg.97",
        "id": "J558.7pg.97*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-10*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-10",
        "length": 319,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TTGTCCACTCCCCAGCCCAGCTGCAGCATTCTGGGGATGACATGGAGGAGCCTGGGTCCTCAGTGAAGTTTTCCTGCATGGCTTCTGGCTATACCTTCACTGACCACTCTATGCATGTGGTAAAACAATGAAACAGAGGCCTGGACAGGGACTGGAGTGGATTGAATGGGTTGGACCTACATGTGGTGGTACTGTATATGCTAGGAAGTTTCAAGGCAAGGCCACACTGACTGTAGACAAATCAGCCATCACAACCTACATGCAGCTCAGTAGCCTGACATCTGAAAACTCTGCAGTCTATTACTGTGCCATGCAAGGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1472910,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1458872,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1458865
            },
            "3 nonamer": {
                "end": 1458904,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCA",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1458895
            },
            "end": 1458865,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGCCTACTCCCAGATCCAGCTGCAACAGTCAGGAGCTGAGCTGGCGAGTCCTGGGGCATCAGTGACACTGTCCTGCAAGGCTTCTGGCTACACATTTACTGACCATATTATGAATTGGGTAAAAAAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAAGGATTTATCCAGTAAGTGGTGAAACTAACTACAATCAAAAGTTCATGGGCAAGGCCACATTCTCTGTAGACCGGTCCTCCAGCACAGTGTACATGGTGTTGAACAGCCTGACATCTGAGGACCCTGCTGTCTATTACTGTGGAAGG",
                "read frame": 1,
                "strand": True
            },
            "start": 1458560
        },
        "accession": "BN000872",
        "allele": "J558.8.98*01",
        "confirmed": True,
        "end": 1458865,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.8.98",
        "id": "J558.8.98*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-11*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-11",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGCCTACTCCCAGATCCAGCTGCAACAGTCAGGAGCTGAGCTGGCGAGTCCTGGGGCATCAGTGACACTGTCCTGCAAGGCTTCTGGCTACACATTTACTGACCATATTATGAATTGGGTAAAAAAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAAGGATTTATCCAGTAAGTGGTGAAACTAACTACAATCAAAAGTTCATGGGCAAGGCCACATTCTCTGTAGACCGGTCCTCCAGCACAGTGTACATGGTGTTGAACAGCCTGACATCTGAGGACCCTGCTGTCTATTACTGTGGAAGG",
            "read frame": 1,
            "strand": True
        },
        "start": 1458560,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 210192,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 210185
            },
            "3 nonamer": {
                "end": 210224,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACAC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 210215
            },
            "end": 210185,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTGCCAGGTCCAGCTGAAGCAGTCTGGAGCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACACCTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAAAGATTGGTCCTGGAAGTGGTAGTACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 209880
        },
        "accession": "BN000872",
        "allele": "J558.80.186*01",
        "confirmed": True,
        "end": 210185,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.80.186",
        "id": "J558.80.186*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-77*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-77",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTGCCAGGTCCAGCTGAAGCAGTCTGGAGCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACACCTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGCCTTGAGTGGATTGGAAAGATTGGTCCTGGAAGTGGTAGTACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 209880,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 203286,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 203279
            },
            "3 nonamer": {
                "end": 203318,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 203309
            },
            "end": 203279,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTTCAGCTGCAACAGTCTGACGCTGAGTTGGTGAAACCTGGAGCTTCAGTGAAGATATCCTGCAAGGTTTCTGGCTACACCTTCACTGACCATACTATTCACTGGATGAAGCAGAGGCCTGAACAGGGCCTGGAATGGATTGGATATATTTATCCTAGAGATGGTAGTACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAACAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 202974
        },
        "accession": "BN000872",
        "allele": "J558.81.187*01",
        "confirmed": True,
        "end": 203279,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.81.187",
        "id": "J558.81.187*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-78*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-78",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTTCAGCTGCAACAGTCTGACGCTGAGTTGGTGAAACCTGGAGCTTCAGTGAAGATATCCTGCAAGGTTTCTGGCTACACCTTCACTGACCATACTATTCACTGGATGAAGCAGAGGCCTGAACAGGGCCTGGAATGGATTGGATATATTTATCCTAGAGATGGTAGTACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAACAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 202974,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 183963,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 183956
            },
            "3 nonamer": {
                "end": 183995,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 183986
            },
            "end": 183956,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCTCAGGTTCAGCTGCAGTAGTCTGGAGCTGAACTGGTGAAACCAGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACAAACTACATCATAAACTGGGTAAAGGAAGGGCCTGGACATGGCCTTGAGTGGATTGGATGGATTTCTCCTGAATATGGTCATACTTACTACAATCAGAAGTTCAAGGGCAAGGCCACATTCACTGCAGACACATCCTCCAGCACAGCCTACATGGAGCTCAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 183651
        },
        "accession": "BN000872",
        "allele": "J558.82pg.188*01",
        "confirmed": True,
        "end": 183956,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.82pg.188",
        "id": "J558.82pg.188*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-79*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-79",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCTCAGGTTCAGCTGCAGTAGTCTGGAGCTGAACTGGTGAAACCAGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACAAACTACATCATAAACTGGGTAAAGGAAGGGCCTGGACATGGCCTTGAGTGGATTGGATGGATTTCTCCTGAATATGGTCATACTTACTACAATCAGAAGTTCAAGGGCAAGGCCACATTCACTGCAGACACATCCTCCAGCACAGCCTACATGGAGCTCAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 183651,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 159960,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 159953
            },
            "3 nonamer": {
                "end": 159992,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAATCC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 159983
            },
            "end": 159953,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAATCCCAGGTTCAGCTGCAGCAGTCTGGGGCTGAGCTGGTGAAGCCTGGGGCCTCAGTGAAGATTTCCTGCAAAGCTTCTGGCTACGCATTCAGTAGCTACTGGATGAACTGGGTGAAGCAGAGGCCTGGAAAGGGTCTTGAGTGGATTGGACAGATTTATCCTGGAGATGGTGATACTAACTACAACGGAAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 159648
        },
        "accession": "BN000872",
        "allele": "J558.83.189*01",
        "confirmed": True,
        "end": 159953,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.83.189",
        "id": "J558.83.189*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-80*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-80",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAATCCCAGGTTCAGCTGCAGCAGTCTGGGGCTGAGCTGGTGAAGCCTGGGGCCTCAGTGAAGATTTCCTGCAAAGCTTCTGGCTACGCATTCAGTAGCTACTGGATGAACTGGGTGAAGCAGAGGCCTGGAAAGGGTCTTGAGTGGATTGGACAGATTTATCCTGGAGATGGTGATACTAACTACAACGGAAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 159648,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 152023,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 152016
            },
            "3 nonamer": {
                "end": 152055,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 152046
            },
            "end": 152016,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCAATCCCAGGTTCAGCTGCAGCAGTCTGGAGCTGAGCTGGCGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACAAGCTATGGTATAAGCTGGGTGAAGCAGAGAACTGGACAGGGCCTTGAGTGGATTGGAGAGATTTATCCTAGAAGTGGTAATACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCGTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 151711
        },
        "accession": "BN000872",
        "allele": "J558.84.190*01",
        "confirmed": True,
        "end": 152016,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.84.190",
        "id": "J558.84.190*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-81*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-81",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAATCCCAGGTTCAGCTGCAGCAGTCTGGAGCTGAGCTGGCGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACAAGCTATGGTATAAGCTGGGTGAAGCAGAGAACTGGACAGGGCCTTGAGTGGATTGGAGAGATTTATCCTAGAAGTGGTAATACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCGTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 151711,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 119773,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 119766
            },
            "3 nonamer": {
                "end": 119804,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 119795
            },
            "end": 119766,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCCTCAGTGAAGATTTCCTGCAAGGCTTCTGGCTACGCATTCAGTAGCTCCTGGATGAACTGGGTGAAGCAGAGGCCTGGAAAGGGTCTTGAGTGGATTGGACGGATTTATCCTGGAGATGGAGATACTAACTACAATGGGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAACTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTACTTCTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 119461
        },
        "accession": "BN000872",
        "allele": "J558.85.191*01",
        "confirmed": True,
        "end": 119766,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.85.191",
        "id": "J558.85.191*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-82*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-82",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCCTCAGTGAAGATTTCCTGCAAGGCTTCTGGCTACGCATTCAGTAGCTCCTGGATGAACTGGGTGAAGCAGAGGCCTGGAAAGGGTCTTGAGTGGATTGGACGGATTTATCCTGGAGATGGAGATACTAACTACAATGGGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACAAATCCTCCAGCACAGCCTACATGCAACTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTACTTCTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 119461,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 108488,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 108481
            },
            "3 nonamer": {
                "end": 108520,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAATCC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 108511
            },
            "end": 108481,
            "length": 291,
            "sequence": {
                "nucleotide": "GTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACATTCACTGACTACTACATGCACTGGGTGAAGCAGAAGCCTGGGAAGGGCCTTGAGTGGATTGGAGAGATTTATCCTGGAAGTGGTAATACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 108190
        },
        "accession": "BN000872",
        "allele": "J558.86.192*01",
        "confirmed": True,
        "end": 108481,
        "family": "J558",
        "framed": True,
        "functionality": "P",
        "gene": "J558.86.192",
        "id": "J558.86.192*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-83*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-83",
        "length": 291,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTCCAGCTGCAGCAGTCTGGACCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGATGTCCTGCAAGGCTTCTGGATACACATTCACTGACTACTACATGCACTGGGTGAAGCAGAAGCCTGGGAAGGGCCTTGAGTGGATTGGAGAGATTTATCCTGGAAGTGGTAATACTTACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGCAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTTCTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 108190,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 91603,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 91596
            },
            "3 nonamer": {
                "end": 91635,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 91626
            },
            "end": 91596,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCATTGCCAGATCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACACCTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTATCCTGGAAGCGGTAATACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 91291
        },
        "accession": "BN000872",
        "allele": "J558.87.193*01",
        "confirmed": True,
        "end": 91596,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.87.193",
        "id": "J558.87.193*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-84*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-84",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCATTGCCAGATCCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGATATCCTGCAAGGCTTCTGGCTACACCTTCACTGACTACTATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTATCCTGGAAGCGGTAATACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 91291,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 72277,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 72270
            },
            "3 nonamer": {
                "end": 72309,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAATCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 72300
            },
            "end": 72270,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGTTGTCCTGCAAGGCTTCTGGCTACACCTTCACAAGCTACGATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTATCCTAGAGATGGTAGTACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCGTACATGGAGCTCCACAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 71965
        },
        "accession": "BN000872",
        "allele": "J558.88.194*01",
        "confirmed": True,
        "end": 72270,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.88.194",
        "id": "J558.88.194*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-85*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-85",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGTTCAGCTGCAGCAGTCTGGACCTGAGCTGGTGAAGCCTGGGGCTTCAGTGAAGTTGTCCTGCAAGGCTTCTGGCTACACCTTCACAAGCTACGATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGACTTGAGTGGATTGGATGGATTTATCCTAGAGATGGTAGTACTAAGTACAATGAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACACATCCTCCAGCACAGCGTACATGGAGCTCCACAGCCTGACATCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 71965,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 62646,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 62639
            },
            "3 nonamer": {
                "end": 62677,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 62668
            },
            "end": 62639,
            "length": 307,
            "sequence": {
                "nucleotide": "GTGTGCACTCCCAGGTCCAGCTGCAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCCTCAGTGAAGATTTCCTGCAAGGCTTTTGGCTACACCTTCACAAACCATCATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGCCTGGACTGGATTGGATATATTAATCCTTATAATGATTATACTAGCTACAGAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTATATGGAGCTTAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 62332
        },
        "accession": "BN000872",
        "allele": "J558.89pg.195*01",
        "confirmed": True,
        "end": 62639,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "J558.89pg.195",
        "id": "J558.89pg.195*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-86*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-86",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTGCACTCCCAGGTCCAGCTGCAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCCTCAGTGAAGATTTCCTGCAAGGCTTTTGGCTACACCTTCACAAACCATCATATAAACTGGGTGAAGCAGAGGCCTGGACAGGGCCTGGACTGGATTGGATATATTAATCCTTATAATGATTATACTAGCTACAGAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTATATGGAGCTTAGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 62332,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1455265,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1455258
            },
            "3 nonamer": {
                "end": 1455297,
                "length": 9,
                "sequence": {
                    "nucleotide": "CAGAAACCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1455288
            },
            "end": 1455258,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGTCCACTCCCAGGCTTATCTACAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCCTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACATTTACCAGTTACAATATGCACTGGGTAAAGCAGACACCTAGACAGGGCCTGGAATGGATTGGAGCTATTTATCCAGGAAATGGTGATACTTCCTACAATCAGAAGTTCAAGGGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAAGACTCTGCGGTCTATTTCTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1454953
        },
        "accession": "BN000872",
        "allele": "J558.9.99*01",
        "confirmed": True,
        "end": 1455258,
        "family": "J558",
        "framed": True,
        "functionality": "F",
        "gene": "J558.9.99",
        "id": "J558.9.99*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-12*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-12",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCACTCCCAGGCTTATCTACAGCAGTCTGGGGCTGAGCTGGTGAGGCCTGGGGCCTCAGTGAAGATGTCCTGCAAGGCTTCTGGCTACACATTTACCAGTTACAATATGCACTGGGTAAAGCAGACACCTAGACAGGGCCTGGAATGGATTGGAGCTATTTATCCAGGAAATGGTGATACTTCCTACAATCAGAAGTTCAAGGGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAAGACTCTGCGGTCTATTTCTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1454953,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 1719406,
            "length": 304,
            "sequence": {
                "nucleotide": "GTATCCATTTCCAAATCCAACTGCAACCGTCTGGGGCCGACCTGATTAATCCTGGATCCTCAATGAAGGTGTCTTGCAAGGTTCCTGGCTACAGTTTCACTAGCTACTATATGACCTGGGTAAAAGAGAGTCCTGGACAGGGTCTAGAATAGATTAGAGAAACCACCCTAGTAGTTGCAGTATAAGTTATGCACAGAAGTTTCAAGGACACATCTCCATGACTAGGCACATATCCTCCAGCATGGCCTACATGGAGCTCAGCAGGCTGACCTCTGAGGACACTTCTGTTTATTGCTGTTTATTG",
                "read frame": 0,
                "strand": True
            },
            "start": 1719102
        },
        "accession": "BN000872",
        "allele": "PG.16.76*01",
        "confirmed": True,
        "end": 1719404,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "PG.16.76",
        "id": "PG.16.76*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-1*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-1",
        "length": 291,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAAATCCAACTGCAACCGTCTGGGGCCGACCTGATTAATCCTGGATCCTCAATGAAGGTGTCTTGCAAGGTTCCTGGCTACAGTTTCACTAGCTACTATATGACCTGGGTAAAAGAGAGTCCTGGACAGGGTCTAGAATAGATTAGAGAAACCACCCTAGTAGTTGCAGTATAAGTTATGCACAGAAGTTTCAAGGACACATCTCCATGACTAGGCACATATCCTCCAGCATGGCCTACATGGAGCTCAGCAGGCTGACCTCTGAGGACACTTCTGTTTATTGCTGTTTAT",
            "read frame": 0,
            "strand": True
        },
        "start": 1719113,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 1589438,
            "length": 247,
            "sequence": {
                "nucleotide": "CAAGGTCCAGCTGCAGCAGTCTGGAACTGAGTTTGGGAGGTCTGGAGCCTCAGAGAAGTAGTCCTGCAAGTCTTCTGGTTACAAATTCACCAACTATTATATGCACCGGATAGGGCAGAATCGTGGAAATAGCCTAGACTAGATTGGATATATTTATCGTGGAATTGGCCATCCTGGCTATGCCCACAACTTCAGGGATAAGACCACACTGACTTCAGACAAATCCTCTAGCACAGCCTACATGGAA",
                "read frame": 2,
                "strand": True
            },
            "start": 1589191
        },
        "accession": "BN000872",
        "allele": "PG.17.87*01",
        "confirmed": True,
        "end": 1589494,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "PG.17.87",
        "id": "PG.17.87*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-3*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-3",
        "length": 302,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "AAGGTCCAGCTGCAGCAGTCTGGAACTGAGTTTGGGAGGTCTGGAGCCTCAGAGAAGTAGTCCTGCAAGTCTTCTGGTTACAAATTCACCAACTATTATATGCACCGGATAGGGCAGAATCGTGGAAATAGCCTAGACTAGATTGGATATATTTATCGTGGAATTGGCCATCCTGGCTATGCCCACAACTTCAGGGATAAGACCACACTGACTTCAGACAAATCCTCTAGCACAGCCTACATGGAATTTCACATCCATTACTCCAGGAATCTTCTTGGGGTATCTCCCCATGACGTAACTCA",
            "read frame": 0,
            "strand": True
        },
        "start": 1589192,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 1538091,
            "length": 244,
            "sequence": {
                "nucleotide": "TCCAGCTGCAGCAGTCTGGAACTGAGTTTGGGTGGACTAGAGCCTCAGGGAAGTAGTCCTGCAAGTCTTCTGGTTACAAATTCACCAACTACTATATGCACTGGATAGGGCAGAATCGTGGAAATAGCCTAGACTAGATTGGATATATTTATCGTGGAATTGGCCATCCTGGCTATGCTCACAACTTCAAGGACAAGACCACACTGACTTCAGACAAATCCTCTAGCACAGCCTACATGGAATT",
                "read frame": 1,
                "strand": True
            },
            "start": 1537847
        },
        "accession": "BN000872",
        "allele": "PG.18.92*01",
        "confirmed": True,
        "end": 1538136,
        "family": "J558",
        "framed": False,
        "functionality": "P",
        "gene": "PG.18.92",
        "id": "PG.18.92*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV1-6*01",
        "imgt family name": "IGHV1",
        "imgt gene name": "IGHV1-6",
        "length": 293,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGGTCCAGCTGCAGCAGTCTGGAACTGAGTTTGGGTGGACTAGAGCCTCAGGGAAGTAGTCCTGCAAGTCTTCTGGTTACAAATTCACCAACTACTATATGCACTGGATAGGGCAGAATCGTGGAAATAGCCTAGACTAGATTGGATATATTTATCGTGGAATTGGCCATCCTGGCTATGCTCACAACTTCAAGGACAAGACCACACTGACTTCAGACAAATCCTCTAGCACAGCCTACATGGAATTTCACATCCATTACTCCAGGAATCTTCTTGGGGTATCTCCCCATGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1537843,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2467782,
            "length": 283,
            "sequence": {
                "nucleotide": "TGTGAGGTGAAGCTAGTGGAGTCTAGGGGAGGCTTAATGCAGCCTGAAAGGTCCGTGATACTCTCCTGTGCTGCCTCTGGATGCACTGTCAGTGACTACTGGTTGGCCTGGGTTTGCCAGGCTCCAAAGAAGGGGCTGGATGGGGGGGGGCATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGAAGAAGTGGGCTGGACATAACATATTCAGAATCAATATTTAGAGATTCTAATCCTTAAAGACATCACTTTGACCAAGTATATATGC",
                "read frame": 1,
                "strand": True
            },
            "start": 2467499
        },
        "accession": "BN000872",
        "allele": "7183.5pg.7*01",
        "confirmed": True,
        "end": 2467796,
        "family": "J606",
        "framed": False,
        "functionality": "P",
        "gene": "7183.5pg.7",
        "id": "7183.5pg.7*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-1*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-1",
        "length": 294,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GAGGTGAAGCTAGTGGAGTCTAGGGGAGGCTTAATGCAGCCTGAAAGGTCCGTGATACTCTCCTGTGCTGCCTCTGGATGCACTGTCAGTGACTACTGGTTGGCCTGGGTTTGCCAGGCTCCAAAGAAGGGGCTGGATGGGGGGGGGCATTAATTTTTCATGGTGGTGGTAGCACCTCCTATGCAGACACCTTGAAGAAGTGGGCTGGACATAACATATTCAGAATCAATATTTAGAGATTCTAATCCTTAAAGACATCACTTTGACCAAGTATATATGCAACATGTTATGGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 2467502,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1679404,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1679397
            },
            "3 nonamer": {
                "end": 1679436,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1679427
            },
            "end": 1679397,
            "length": 310,
            "sequence": {
                "nucleotide": "GGGTCCAGAGTGAAGTGAAGCTTGAGGAGTCTGGAGGAGGCTTGGTGCAACCTGGAGGATCCATGAAACTCTCCTGTGTTGCCTCTGGATTCACTTTCAGTAACTACTGGATGAACTGGGTCCGCCAGTCTCCAGAGAAGGGGCTTGAGTGGGTTGCTCAAATTAGATTGAAATCTGATAATTATGCAACACATTATGCGGAGTCTGTGAAAGGGAGGTTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTCTACCTGCAAATGAACAACTTAAGGGCTGAAGACACTGGAATTTATTACTGCACAGG",
                "read frame": 1,
                "strand": True
            },
            "start": 1679087
        },
        "accession": "BN000872",
        "allele": "J606.1.79*01",
        "confirmed": True,
        "end": 1679397,
        "family": "J606",
        "framed": True,
        "functionality": "F",
        "gene": "J606.1.79",
        "id": "J606.1.79*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-3*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-3",
        "length": 310,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGGTCCAGAGTGAAGTGAAGCTTGAGGAGTCTGGAGGAGGCTTGGTGCAACCTGGAGGATCCATGAAACTCTCCTGTGTTGCCTCTGGATTCACTTTCAGTAACTACTGGATGAACTGGGTCCGCCAGTCTCCAGAGAAGGGGCTTGAGTGGGTTGCTCAAATTAGATTGAAATCTGATAATTATGCAACACATTATGCGGAGTCTGTGAAAGGGAGGTTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTCTACCTGCAAATGAACAACTTAAGGGCTGAAGACACTGGAATTTATTACTGCACAGG",
            "read frame": 1,
            "strand": True
        },
        "start": 1679087,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1664642,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTA",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1664635
            },
            "3 nonamer": {
                "end": 1664674,
                "length": 9,
                "sequence": {
                    "nucleotide": "CCACCAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1664665
            },
            "end": 1664635,
            "length": 311,
            "sequence": {
                "nucleotide": "GTGTCCAGAGTGATGTGAACCTGGAAGTGTCTGGAGGAGGCTTAGTTAAACCTGGAGGATCCATGCAACTCTTTTGTGTAGCCTCTGGATTTACTTTTGTAGATGGCTGGATGGACTGGGTCCGCCAGTCTCCAGAGAAGGGTCTTGAGTGGGTTGCTGAAATTGCAAACAAAGCTAATAATTATGCAACATATTATCCCGAGTCTGTGAAAGGCAGATTCACCATCTCAAGAGATGATTTCAAAAGTAGTGTCTACCTGCATATGAACAGCTTAAGAGCCGAAGATACAGGCATTTATTACTGTACAAGG",
                "read frame": 1,
                "strand": True
            },
            "start": 1664324
        },
        "accession": "BN000872",
        "allele": "J606.2.80*01",
        "confirmed": True,
        "end": 1664635,
        "family": "J606",
        "framed": True,
        "functionality": "F",
        "gene": "J606.2.80",
        "id": "J606.2.80*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-4*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-4",
        "length": 311,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGAGTGATGTGAACCTGGAAGTGTCTGGAGGAGGCTTAGTTAAACCTGGAGGATCCATGCAACTCTTTTGTGTAGCCTCTGGATTTACTTTTGTAGATGGCTGGATGGACTGGGTCCGCCAGTCTCCAGAGAAGGGTCTTGAGTGGGTTGCTGAAATTGCAAACAAAGCTAATAATTATGCAACATATTATCCCGAGTCTGTGAAAGGCAGATTCACCATCTCAAGAGATGATTTCAAAAGTAGTGTCTACCTGCATATGAACAGCTTAAGAGCCGAAGATACAGGCATTTATTACTGTACAAGG",
            "read frame": 1,
            "strand": True
        },
        "start": 1664324,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1654520,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1654513
            },
            "3 nonamer": {
                "end": 1654552,
                "length": 9,
                "sequence": {
                    "nucleotide": "CCACAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1654543
            },
            "end": 1654513,
            "length": 311,
            "sequence": {
                "nucleotide": "GTGTCCAGAGTGAAGTGAAAATTGAGGAGTCAGGAGGAGGCTTGGTCCAACCTGGAGGATCCATGAAACTCTCTTGTGCAGCCTCTGGATTCACTTTCAGTGATTACAGGATGGACTGGGTCCACCACTCTACAGAGAATGGGTTGGAGTGGGTTGCTGAAATTAGAAACAAAGCTAGTAATTATGCAACATATTATGTGGAGTCTGTGAATGGGAGGTTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTCTACCTGCAAATGAACAGCTTAAGAGCTGAAGATACTGGCATTTATTACTGTACAAGG",
                "read frame": 1,
                "strand": True
            },
            "start": 1654202
        },
        "accession": "BN000872",
        "allele": "J606.3.81*01",
        "confirmed": True,
        "end": 1654513,
        "family": "J606",
        "framed": True,
        "functionality": "F",
        "gene": "J606.3.81",
        "id": "J606.3.81*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-5*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-5",
        "length": 311,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGAGTGAAGTGAAAATTGAGGAGTCAGGAGGAGGCTTGGTCCAACCTGGAGGATCCATGAAACTCTCTTGTGCAGCCTCTGGATTCACTTTCAGTGATTACAGGATGGACTGGGTCCACCACTCTACAGAGAATGGGTTGGAGTGGGTTGCTGAAATTAGAAACAAAGCTAGTAATTATGCAACATATTATGTGGAGTCTGTGAATGGGAGGTTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTCTACCTGCAAATGAACAGCTTAAGAGCTGAAGATACTGGCATTTATTACTGTACAAGG",
            "read frame": 1,
            "strand": True
        },
        "start": 1654202,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1636326,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1636319
            },
            "3 nonamer": {
                "end": 1636358,
                "length": 9,
                "sequence": {
                    "nucleotide": "CCACAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1636349
            },
            "end": 1636319,
            "length": 311,
            "sequence": {
                "nucleotide": "GTGTCCAGAGTGAAGTGAAGCTTGAGGAGTCTGGAGGAGGCTTGGTGCAACCTGGAGGATCCATGAAACTCTCTTGTGCTGCCTCTGGATTCACTTTTAGTGACGCCTGGATGGACTGGGTCCGCCAGTCTCCAGAGAAGGGGCTTGAGTGGGTTGCTGAAATTAGAAACAAAGCTAATAATCATGCAACATACTATGCTGAGTCTGTGAAAGGGAGGTTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTCTACCTGCAAATGAACAGCTTAAGAGCTGAAGACACTGGCATTTATTACTGTACCAGG",
                "read frame": 0,
                "strand": True
            },
            "start": 1636008
        },
        "accession": "BN000872",
        "allele": "J606.4.82*01",
        "confirmed": True,
        "end": 1636319,
        "family": "J606",
        "framed": True,
        "functionality": "F",
        "gene": "J606.4.82",
        "id": "J606.4.82*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-6*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-6",
        "length": 311,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGAGTGAAGTGAAGCTTGAGGAGTCTGGAGGAGGCTTGGTGCAACCTGGAGGATCCATGAAACTCTCTTGTGCTGCCTCTGGATTCACTTTTAGTGACGCCTGGATGGACTGGGTCCGCCAGTCTCCAGAGAAGGGGCTTGAGTGGGTTGCTGAAATTAGAAACAAAGCTAATAATCATGCAACATACTATGCTGAGTCTGTGAAAGGGAGGTTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTCTACCTGCAAATGAACAGCTTAAGAGCTGAAGACACTGGCATTTATTACTGTACCAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1636008,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1615490,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1615483
            },
            "3 nonamer": {
                "end": 1615522,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1615513
            },
            "end": 1615483,
            "length": 311,
            "sequence": {
                "nucleotide": "GTGTCCAGTGTGAGGAGAAGCTGGATGAGTCTGGAGGAGGCTTGGTGCAACCTGGGAGGTCCATGAAACTCTCCTGTGTTGCCTCTGGATTCACTTTTACTAACTCCTGGATGAACTGGTTCTGCCAGTCTCCAGAGAAAGGACTGGAGTGGGTAGCACAAATTAAAAGCAAACCTTATAATTATGAAACATATTATTCAGATTCTGTGAAAGGCAGATTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTATACCTGCAAATGAACAACTTAAGAGCTGAAGACACGGGCATCTATTACTGTACATGG",
                "read frame": 1,
                "strand": True
            },
            "start": 1615172
        },
        "accession": "BN000872",
        "allele": "J606.5.83*01",
        "confirmed": True,
        "end": 1615483,
        "family": "J606",
        "framed": True,
        "functionality": "F",
        "gene": "J606.5.83",
        "id": "J606.5.83*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV6-7*01",
        "imgt family name": "IGHV6",
        "imgt gene name": "IGHV6-7",
        "length": 311,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCAGTGTGAGGAGAAGCTGGATGAGTCTGGAGGAGGCTTGGTGCAACCTGGGAGGTCCATGAAACTCTCCTGTGTTGCCTCTGGATTCACTTTTACTAACTCCTGGATGAACTGGTTCTGCCAGTCTCCAGAGAAAGGACTGGAGTGGGTAGCACAAATTAAAAGCAAACCTTATAATTATGAAACATATTATTCAGATTCTGTGAAAGGCAGATTCACCATCTCAAGAGATGATTCCAAAAGTAGTGTATACCTGCAAATGAACAACTTAAGAGCTGAAGACACGGGCATCTATTACTGTACATGG",
            "read frame": 1,
            "strand": True
        },
        "start": 1615172,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 1980945,
            "length": 288,
            "sequence": {
                "nucleotide": "TGCAAAACTCCATATTTCTCTCATTCCTTCTTCTCCATCAAATATTTGTTCATCATAGCTGTTTCTGCTACCTACTTGGCAAGTGGATGAAGACATAGAACTCTCTGCTCCTTCTGTGCCATGCCTGCCTAGATACTACCATGCTCCCCTCTTGATGATAATGGACTTAACCTCTGAACCTGTAAGCCAGCCCCAATTAAATGTTGTTTTTTATAAGACTTGCTTTGGTGATGGTGTCTGTTCACAGCAGTAAAACCCTAACTAAGACAGCAAGCTAGGTTGAGTAGC",
                "read frame": 0,
                "strand": True
            },
            "start": 1980657
        },
        "accession": "BN000872",
        "allele": "PG.10.56*01",
        "confirmed": True,
        "end": 1980945,
        "family": "J606",
        "framed": True,
        "functionality": "P",
        "gene": "PG.10.56",
        "id": "PG.10.56*01-C57BL/6",
        "identified": True,
        "length": 288,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGCAAAACTCCATATTTCTCTCATTCCTTCTTCTCCATCAAATATTTGTTCATCATAGCTGTTTCTGCTACCTACTTGGCAAGTGGATGAAGACATAGAACTCTCTGCTCCTTCTGTGCCATGCCTGCCTAGATACTACCATGCTCCCCTCTTGATGATAATGGACTTAACCTCTGAACCTGTAAGCCAGCCCCAATTAAATGTTGTTTTTTATAAGACTTGCTTTGGTGATGGTGTCTGTTCACAGCAGTAAAACCCTAACTAAGACAGCAAGCTAGGTTGAGTAGC",
            "read frame": 0,
            "strand": True
        },
        "start": 1980657,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2275018,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2275011
            },
            "3 nonamer": {
                "end": 2275050,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAATACT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2275041
            },
            "end": 2275011,
            "length": 304,
            "sequence": {
                "nucleotide": "GTATCCTTTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACTGTCTCTGGGTTCTCATTAACCAGCTATGGTGTAGACTGGGTTCGCCAGTCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTAATATGGGGTGTTGGAAGCACAAATTATAATTCAGCTCTCAAATCCAGACTGAGCATCAGCAAGGACAACTCCAAGAGCCAAGTTTTCTTAAAAATGAACAGTCTGCAAACTGATGACACAGCCATGTACTACTGTGCCAGTGA",
                "read frame": 1,
                "strand": True
            },
            "start": 2274707
        },
        "accession": "BN000872",
        "allele": "Q52.10.33*01",
        "confirmed": True,
        "end": 2275011,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.10.33",
        "id": "Q52.10.33*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-6-8*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-6-8",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCTTTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACTGTCTCTGGGTTCTCATTAACCAGCTATGGTGTAGACTGGGTTCGCCAGTCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTAATATGGGGTGTTGGAAGCACAAATTATAATTCAGCTCTCAAATCCAGACTGAGCATCAGCAAGGACAACTCCAAGAGCCAAGTTTTCTTAAAAATGAACAGTCTGCAAACTGATGACACAGCCATGTACTACTGTGCCAGTGA",
            "read frame": 1,
            "strand": True
        },
        "start": 2274707,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2263842,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2263835
            },
            "3 nonamer": {
                "end": 2263874,
                "length": 9,
                "sequence": {
                    "nucleotide": "AACAAAAGC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2263865
            },
            "end": 2263835,
            "length": 303,
            "sequence": {
                "nucleotide": "GTGTCCTGTCACAAGTGCAGATGAAGGAGTCAGGACCTGACCTTGTGCAGCCATCACAGACTCTGTCTCTCACCTGCACTGTCTCTGGGTTCTCATTAAGTAGCTATGGTGTACATTGGTTTCGCAAGCCTCCGAGAAAGGGATTGGAATGGTTGGGAGGAATATGGTCTGGTGGAAGCATATACTATACTCCAGCTCTCAGTTCCCGACTGAGTGTCAGCAGGGACACCTCTAAGAGCCAAGTTTTCTTTAAAATGAGCAGTCTGCAAAGTGAAGACACGGCTGTGTACCACTGTGCCAGATA",
                "read frame": 2,
                "strand": True
            },
            "start": 2263531
        },
        "accession": "BN000872",
        "allele": "Q52.11.34*01",
        "confirmed": True,
        "end": 2263835,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.11.34",
        "id": "Q52.11.34*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-7*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-7",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTGTCACAAGTGCAGATGAAGGAGTCAGGACCTGACCTTGTGCAGCCATCACAGACTCTGTCTCTCACCTGCACTGTCTCTGGGTTCTCATTAAGTAGCTATGGTGTACATTGGTTTCGCAAGCCTCCGAGAAAGGGATTGGAATGGTTGGGAGGAATATGGTCTGGTGGAAGCATATACTATACTCCAGCTCTCAGTTCCCGACTGAGTGTCAGCAGGGACACCTCTAAGAGCCAAGTTTTCTTTAAAATGAGCAGTCTGCAAAGTGAAGACACGGCTGTGTACCACTGTGCCAGATA",
            "read frame": 2,
            "strand": True
        },
        "start": 2263531,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2193794,
            "length": 297,
            "sequence": {
                "nucleotide": "GAGGTACAACTGAAGTAGTCATGACCTGACCTGCTGCAATCCTCACAGGCTGTCTCTCATCTGCACTGTCTCTGGCTTCTCATTAACCATCTATGGTGTAAACTGGGTCCTCCAGCCACTAGGAAGGGGATTGCAGAGGATGCCAGCAATATGGAGAGGTGAAAGCACAGAGTATAATTCAGATCTCAAATCCTGAATCAGCATCAGTAGGGACACATCCAAGAGTCAACTTTTCTTAAAACTGAACAGAATTCAAACTGAGGACATAGCCATCTATGACCCTGTCAGAGAAACACT",
                "read frame": 1,
                "strand": True
            },
            "start": 2193497
        },
        "accession": "BN000872",
        "allele": "Q52.12pg.39*01",
        "confirmed": True,
        "end": 2193788,
        "family": "Q52",
        "framed": False,
        "functionality": "P",
        "gene": "Q52.12pg.39",
        "id": "Q52.12pg.39*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-8*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-8",
        "length": 299,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ACCCTGTAGAGGTACAACTGAAGTAGTCATGACCTGACCTGCTGCAATCCTCACAGGCTGTCTCTCATCTGCACTGTCTCTGGCTTCTCATTAACCATCTATGGTGTAAACTGGGTCCTCCAGCCACTAGGAAGGGGATTGCAGAGGATGCCAGCAATATGGAGAGGTGAAAGCACAGAGTATAATTCAGATCTCAAATCCTGAATCAGCATCAGTAGGGACACATCCAAGAGTCAACTTTTCTTAAAACTGAACAGAATTCAAACTGAGGACATAGCCATCTATGACCCTGTCAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2193489,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2192060,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2192053
            },
            "3 nonamer": {
                "end": 2192092,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2192083
            },
            "end": 2192053,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGTCCTGTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACTTGCACTGTCTCTGGGTTTTCATTAACCAGCTATGGTGTAGACTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTAATATGGGGTGGTGGAAGCACAAATTATAATTCAGCTCTCATGTCCAGACTGAGCATCAGCAAAGACAACTCCAAGAGCCAAGTTTTCTTAAAAATGAACAGTCTGCAAACTGATGACACAGCCATGTACTACTGTGCCAAACA",
                "read frame": 0,
                "strand": True
            },
            "start": 2191749
        },
        "accession": "BN000872",
        "allele": "Q52.13.40*01",
        "confirmed": True,
        "end": 2192053,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.13.40",
        "id": "Q52.13.40*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-9*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-9",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTGTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACTTGCACTGTCTCTGGGTTTTCATTAACCAGCTATGGTGTAGACTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTAATATGGGGTGGTGGAAGCACAAATTATAATTCAGCTCTCATGTCCAGACTGAGCATCAGCAAAGACAACTCCAAGAGCCAAGTTTTCTTAAAAATGAACAGTCTGCAAACTGATGACACAGCCATGTACTACTGTGCCAAACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2191749,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2497015,
            "length": 311,
            "sequence": {
                "nucleotide": "GTGTCCTGTCCCAGTTCCAGCTGAAGCAGTCAGGACCTGCCCTTGTGCCGCCTTCACGGAGCCTGTCCATCACATGCACTGTCTCTGGGTTCTCATTAACCAGCTGTGGTATAAGCTGGCTTCGATGGTCTCCAGGAAACGGAGATAGTGGAATTTTACAACTGTCCCATAGGCCTAGAAAACCCAGTCATCACCCCTGTTCCTCTTCTGCATCCTTACCAATACGTGATCCCCAATCCTTACCCACTATTCATCAAGTATATTATATTTTGCTACTGTTCCCCATTTATACAGGACTCGCCAGATCACAC",
                "read frame": 1,
                "strand": True
            },
            "start": 2496704
        },
        "accession": "BN000872",
        "allele": "Q52.1pg.2*01",
        "confirmed": True,
        "end": 2497009,
        "family": "Q52",
        "framed": False,
        "functionality": "P",
        "gene": "Q52.1pg.2",
        "id": "Q52.1pg.2*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-1*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-1",
        "length": 294,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGTTCCAGCTGAAGCAGTCAGGACCTGCCCTTGTGCCGCCTTCACGGAGCCTGTCCATCACATGCACTGTCTCTGGGTTCTCATTAACCAGCTGTGGTATAAGCTGGCTTCGATGGTCTCCAGGAAACGGAGATAGTGGAATTTTACAACTGTCCCATAGGCCTAGAAAACCCAGTCATCACCCCTGTTCCTCTTCTGCATCCTTACCAATACGTGATCCCCAATCCTTACCCACTATTCATCAAGTATATTATATTTTGCTACTGTTCCCCATTTATACAGGACTCGCCAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2496715,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2482996,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2482989
            },
            "3 nonamer": {
                "end": 2483027,
                "length": 9,
                "sequence": {
                    "nucleotide": "AACAAAAAT",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2483018
            },
            "end": 2482989,
            "length": 302,
            "sequence": {
                "nucleotide": "GTGTCCTATCCCAGGTGCAGCTGAAGCAGTCAGGACCTGGCCTAGTGCAGCCCTCACAGAGCCTGTCCATCACCTGCACAGTCTCTGGTTTCTCATTAACTAGCTATGGTGTACACTGGGTTCGCCAGTCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTGATATGGAGTGGTGGAAGCACAGACTATAATGCAGCTTTCATATCCAGACTGAGCATCAGCAAGGACAATTCCAAGAGCCAAGTTTTCTTTAAAATGAACAGTCTGCAAGCTGATGACACAGCCATATATTACTGTGCCAGAAA",
                "read frame": 1,
                "strand": True
            },
            "start": 2482685
        },
        "accession": "BN000872",
        "allele": "Q52.2.4*01",
        "confirmed": True,
        "end": 2482989,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.2.4",
        "id": "Q52.2.4*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-2*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-2",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTATCCCAGGTGCAGCTGAAGCAGTCAGGACCTGGCCTAGTGCAGCCCTCACAGAGCCTGTCCATCACCTGCACAGTCTCTGGTTTCTCATTAACTAGCTATGGTGTACACTGGGTTCGCCAGTCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTGATATGGAGTGGTGGAAGCACAGACTATAATGCAGCTTTCATATCCAGACTGAGCATCAGCAAGGACAATTCCAAGAGCCAAGTTTTCTTTAAAATGAACAGTCTGCAAGCTGATGACACAGCCATATATTACTGTGCCAGAAA",
            "read frame": 1,
            "strand": True
        },
        "start": 2482685,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2460079,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAATG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2460072
            },
            "3 nonamer": {
                "end": 2460111,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAATACT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2460102
            },
            "end": 2460072,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGTCCTGTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACTGTCTCAGGGTTCTCATTAACCAGCTATGGTGTAAGCTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTAATATGGGGTGACGGGAGCACAAATTATCATTCAGCTCTCATATCCAGACTGAGCATCAGCAAGGATAACTCCAAGAGCCAAGTTTTCTTAAAACTGAACAGTCTGCAAACTGATGACACAGCCACGTACTACTGTGCCAAACC",
                "read frame": 1,
                "strand": True
            },
            "start": 2459768
        },
        "accession": "BN000872",
        "allele": "Q52.3.8*01",
        "confirmed": True,
        "end": 2460072,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.3.8",
        "id": "Q52.3.8*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-3*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-3",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTGTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACTGTCTCAGGGTTCTCATTAACCAGCTATGGTGTAAGCTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTAATATGGGGTGACGGGAGCACAAATTATCATTCAGCTCTCATATCCAGACTGAGCATCAGCAAGGATAACTCCAAGAGCCAAGTTTTCTTAAAACTGAACAGTCTGCAAACTGATGACACAGCCACGTACTACTGTGCCAAACC",
            "read frame": 1,
            "strand": True
        },
        "start": 2459768,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2423778,
            "length": 289,
            "sequence": {
                "nucleotide": "TACTTGCATATAAAAGAGTGAGGACTTGAACTGGTGTAGCCTTCACAGACCCTGCCCCATACCTGCACTGTCTCTGGCTTCTCATTATTCCAGCTACCATGTGCACTGAGTGTGTCAGTCTACAGCAAAGGGTCCACAGTGGATGGGAGCAATGTGTAGTGGGGAAACACAGCATACGGTTCAGCTCTCAACTCCCAACTCGTCATCAATAGGGACACATCTATGGGCAAGTGTTCTTAAAACTGTACAATCTTCAAAAAGGGAAACAGTGATGTACTACTGTGCCAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 2423489
        },
        "accession": "BN000872",
        "allele": "Q52.4pg.12*01",
        "confirmed": True,
        "end": 2423778,
        "family": "Q52",
        "framed": True,
        "functionality": "P",
        "gene": "Q52.4pg.12",
        "id": "Q52.4pg.12*01-C57BL/6",
        "identified": True,
        "length": 288,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ACTTGCATATAAAAGAGTGAGGACTTGAACTGGTGTAGCCTTCACAGACCCTGCCCCATACCTGCACTGTCTCTGGCTTCTCATTATTCCAGCTACCATGTGCACTGAGTGTGTCAGTCTACAGCAAAGGGTCCACAGTGGATGGGAGCAATGTGTAGTGGGGAAACACAGCATACGGTTCAGCTCTCAACTCCCAACTCGTCATCAATAGGGACACATCTATGGGCAAGTGTTCTTAAAACTGTACAATCTTCAAAAAGGGAAACAGTGATGTACTACTGTGCCAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2423490,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2417972,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2417965
            },
            "3 nonamer": {
                "end": 2418004,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAAAATT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2417995
            },
            "end": 2417965,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGTCCTATCCCAGGTGCAGCTGAAGCAGTCAGGACCTGGCCTAGTGCAGCCCTCACAGAGCCTGTCCATCACCTGCACAGTCTCTGGTTTCTCATTAACTAGCTATGGTGTACACTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTGATATGGAGTGGTGGAAGCACAGACTATAATGCTGCTTTCATATCCAGACTGAGCATCAGCAAGGACAACTCCAAGAGCCAAGTTTTCTTTAAAATGAACAGTCTGCAAGCTGATGACACTGCCATATACTACTGTGCCAAAAA",
                "read frame": 0,
                "strand": True
            },
            "start": 2417661
        },
        "accession": "BN000872",
        "allele": "Q52.5.13*01",
        "confirmed": True,
        "end": 2417965,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.5.13",
        "id": "Q52.5.13*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-4*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-4",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTATCCCAGGTGCAGCTGAAGCAGTCAGGACCTGGCCTAGTGCAGCCCTCACAGAGCCTGTCCATCACCTGCACAGTCTCTGGTTTCTCATTAACTAGCTATGGTGTACACTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTGATATGGAGTGGTGGAAGCACAGACTATAATGCTGCTTTCATATCCAGACTGAGCATCAGCAAGGACAACTCCAAGAGCCAAGTTTTCTTTAAAATGAACAGTCTGCAAGCTGATGACACTGCCATATACTACTGTGCCAAAAA",
            "read frame": 0,
            "strand": True
        },
        "start": 2417661,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 2388792,
            "length": 301,
            "sequence": {
                "nucleotide": "AATCGTCTCCTACTTGCATCTGAAAGAGCGAGGACTTGAACTGGTGTAGCCTTCACAGACCCTGCCCCATACCTGCACTGTCTCTGGCTTCTCATTATTCAAGCTACCATGTGCACTGAGTGTGTCAGTCTACAGCAAAGGGTCCACAGTGGATGGGAGCAATGTGTAGTGGGGAAACACAGCATACGGTTCAGCTCTCAACTCCCAACTCGTCATCAATAGGGACACATCTATGGGCAAGTGTTCTTAAAACTGAACAATCTTCAAAAAGGAAAACAGTGATGTACTACTGTGCCAGAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 2388491
        },
        "accession": "BN000872",
        "allele": "Q52.6pg.17*01",
        "confirmed": True,
        "end": 2388792,
        "family": "Q52",
        "framed": True,
        "functionality": "P",
        "gene": "Q52.6pg.17",
        "id": "Q52.6pg.17*01-C57BL/6",
        "identified": True,
        "length": 300,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATCGTCTCCTACTTGCATCTGAAAGAGCGAGGACTTGAACTGGTGTAGCCTTCACAGACCCTGCCCCATACCTGCACTGTCTCTGGCTTCTCATTATTCAAGCTACCATGTGCACTGAGTGTGTCAGTCTACAGCAAAGGGTCCACAGTGGATGGGAGCAATGTGTAGTGGGGAAACACAGCATACGGTTCAGCTCTCAACTCCCAACTCGTCATCAATAGGGACACATCTATGGGCAAGTGTTCTTAAAACTGAACAATCTTCAAAAAGGAAAACAGTGATGTACTACTGTGCCAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 2388492,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2385781,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2385774
            },
            "3 nonamer": {
                "end": 2385813,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAAAATT",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2385804
            },
            "end": 2385774,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGTCCTATCCCAGGTGCAGCTGAAGCAGTCAGGACCTGGCCTAGTGCAGCCCTCACAGAGCCTGTCCATAACCTGCACAGTCTCTGGTTTCTCATTAACTAGCTATGGTGTACACTGGGTTCGCCAGTCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTGATATGGAGAGGTGGAAGCACAGACTACAATGCAGCTTTCATGTCCAGACTGAGCATCACCAAGGACAACTCCAAGAGCCAAGTTTTCTTTAAAATGAACAGTCTGCAAGCTGATGACACTGCCATATACTACTGTGCCAAAAA",
                "read frame": 1,
                "strand": True
            },
            "start": 2385470
        },
        "accession": "BN000872",
        "allele": "Q52.7.18*01",
        "confirmed": True,
        "end": 2385774,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.7.18",
        "id": "Q52.7.18*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-5*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-5",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTATCCCAGGTGCAGCTGAAGCAGTCAGGACCTGGCCTAGTGCAGCCCTCACAGAGCCTGTCCATAACCTGCACAGTCTCTGGTTTCTCATTAACTAGCTATGGTGTACACTGGGTTCGCCAGTCTCCAGGAAAGGGTCTGGAGTGGCTGGGAGTGATATGGAGAGGTGGAAGCACAGACTACAATGCAGCTTTCATGTCCAGACTGAGCATCACCAAGGACAACTCCAAGAGCCAAGTTTTCTTTAAAATGAACAGTCTGCAAGCTGATGACACTGCCATATACTACTGTGCCAAAAA",
            "read frame": 1,
            "strand": True
        },
        "start": 2385470,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2354489,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2354482
            },
            "3 nonamer": {
                "end": 2354521,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAATACT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2354512
            },
            "end": 2354482,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGTCCTGTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACCGTCTCAGGGTTCTCATTAACCAGCTATGGTGTACACTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGTAGTGATATGGAGTGATGGAAGCACAACCTATAATTCAGCTCTCAAATCCAGACTGAGCATCAGCAAGGACAACTCCAAGAGCCAAGTTTTCTTAAAAATGAACAGTCTCCAAACTGATGACACAGCCATGTACTACTGTGCCAGACA",
                "read frame": 0,
                "strand": True
            },
            "start": 2354178
        },
        "accession": "BN000872",
        "allele": "Q52.8.22*01",
        "confirmed": True,
        "end": 2354482,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.8.22",
        "id": "Q52.8.22*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-6*03",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-6",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTCCTGTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACCGTCTCAGGGTTCTCATTAACCAGCTATGGTGTACACTGGGTTCGCCAGCCTCCAGGAAAGGGTCTGGAGTGGCTGGTAGTGATATGGAGTGATGGAAGCACAACCTATAATTCAGCTCTCAAATCCAGACTGAGCATCAGCAAGGACAACTCCAAGAGCCAAGTTTTCTTAAAAATGAACAGTCTCCAAACTGATGACACAGCCATGTACTACTGTGCCAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 2354178,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2301305,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2301298
            },
            "3 nonamer": {
                "end": 2301337,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAATACT",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2301328
            },
            "end": 2301298,
            "length": 304,
            "sequence": {
                "nucleotide": "GTGCCCTGTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACTGTCTCTGGGTTCTCATTAACCAGCTATGCTATAAGCTGGGTTCGCCAGCCACCAGGAAAGGGTCTGGAGTGGCTTGGAGTAATATGGACTGGTGGAGGCACAAATTATAATTCAGCTCTCAAATCCAGACTGAGCATCAGCAAAGACAACTCCAAGAGTCAAGTTTTCTTAAAAATGAACAGTCTGCAAACTGATGACACAGCCAGGTACTACTGTGCCAGAAA",
                "read frame": 0,
                "strand": True
            },
            "start": 2300994
        },
        "accession": "BN000872",
        "allele": "Q52.9.29*01",
        "confirmed": True,
        "end": 2301298,
        "family": "Q52",
        "framed": True,
        "functionality": "F",
        "gene": "Q52.9.29",
        "id": "Q52.9.29*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV2-9-1*01",
        "imgt family name": "IGHV2",
        "imgt gene name": "IGHV2-9-1",
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGCCCTGTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAGCCTGTCCATCACATGCACTGTCTCTGGGTTCTCATTAACCAGCTATGCTATAAGCTGGGTTCGCCAGCCACCAGGAAAGGGTCTGGAGTGGCTTGGAGTAATATGGACTGGTGGAGGCACAAATTATAATTCAGCTCTCAAATCCAGACTGAGCATCAGCAAAGACAACTCCAAGAGTCAAGTTTTCTTAAAAATGAACAGTCTGCAAACTGATGACACAGCCAGGTACTACTGTGCCAGAAA",
            "read frame": 0,
            "strand": True
        },
        "start": 2300994,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2174748,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2174741
            },
            "3 nonamer": {
                "end": 2174780,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2174771
            },
            "end": 2174741,
            "length": 315,
            "sequence": {
                "nucleotide": "GTATCCAGTGTGAGGTGAAGCTGGTGGAATCTGGAGGAGGCTTGGTACAGTCTGGGCGTTCTCTGAGACTCTCCTGTGCAACTTCTGGGTTCACCTTCAGTGATTTCTACATGGAGTGGGTCCGCCAAGCTCCAGGGAAGGGACTGGAGTGGATTGCTGCAAGTAGAAACAAAGCTAATGATTATACAACAGAGTACAGTGCATCTGTGAAGGGTCGGTTCATCGTCTCCAGAGACACTTCCCAAAGCATCCTCTACCTTCAGATGAATGCCCTGAGAGCTGAGGACACTGCCATTTATTACTGTGCAAGAGATGCA",
                "read frame": 0,
                "strand": True
            },
            "start": 2174424
        },
        "accession": "BN000872",
        "allele": "S107.1.42*01",
        "confirmed": True,
        "end": 2174741,
        "family": "S107",
        "framed": True,
        "functionality": "P",
        "gene": "S107.1.42",
        "id": "S107.1.42*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV7-1*01",
        "imgt family name": "IGHV7",
        "imgt gene name": "IGHV7-1",
        "length": 317,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCAGTGTGAGGTGAAGCTGGTGGAATCTGGAGGAGGCTTGGTACAGTCTGGGCGTTCTCTGAGACTCTCCTGTGCAACTTCTGGGTTCACCTTCAGTGATTTCTACATGGAGTGGGTCCGCCAAGCTCCAGGGAAGGGACTGGAGTGGATTGCTGCAAGTAGAAACAAAGCTAATGATTATACAACAGAGTACAGTGCATCTGTGAAGGGTCGGTTCATCGTCTCCAGAGACACTTCCCAAAGCATCCTCTACCTTCAGATGAATGCCCTGAGAGCTGAGGACACTGCCATTTATTACTGTGCAAGAGATGCA",
            "read frame": 0,
            "strand": True
        },
        "start": 2174424,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2159131,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2159124
            },
            "3 nonamer": {
                "end": 2159155,
                "length": 9,
                "sequence": {
                    "nucleotide": "AAACAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2159146
            },
            "end": 2159124,
            "length": 309,
            "sequence": {
                "nucleotide": "GTATCCAGTGTGAGGTGAAGCTGGTGGAGTCTGAAGGAGGCTTGGTACAGCCTGGGGGTTCTCTGAGACTCTCCTGTGCAACTTCTGGGTTCACCTTCACTGATTTCTACATGAACTGGGTCTGCCAGCCTCCAAGGAAGGCACTTGAGTGGTTGGGTTTTATTAGAAACAAAGCTAATGGTTACACAACAGACTACAGTGCATCTATGAAGGGTCGGTTCACCATCTCCAGAGATAATTCCCAAAGCATCCTCTATCTTCAAATGAACACACTGAGAACTGAGGACAGTGCCACTTATTACTGTGCAAGAGATACA",
                "read frame": 2,
                "strand": True
            },
            "start": 2158807
        },
        "accession": "BN000872",
        "allele": "S107.2pg.43*01",
        "confirmed": True,
        "end": 2159124,
        "family": "S107",
        "framed": True,
        "functionality": "F",
        "gene": "S107.2pg.43",
        "id": "S107.2pg.43*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV7-2*01",
        "imgt family name": "IGHV7",
        "imgt gene name": "IGHV7-2",
        "length": 317,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCAGTGTGAGGTGAAGCTGGTGGAGTCTGAAGGAGGCTTGGTACAGCCTGGGGGTTCTCTGAGACTCTCCTGTGCAACTTCTGGGTTCACCTTCACTGATTTCTACATGAACTGGGTCTGCCAGCCTCCAAGGAAGGCACTTGAGTGGTTGGGTTTTATTAGAAACAAAGCTAATGGTTACACAACAGACTACAGTGCATCTATGAAGGGTCGGTTCACCATCTCCAGAGATAATTCCCAAAGCATCCTCTATCTTCAAATGAACACACTGAGAACTGAGGACAGTGCCACTTATTACTGTGCAAGAGATACA",
            "read frame": 2,
            "strand": True
        },
        "start": 2158807,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1917939,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1917932
            },
            "3 nonamer": {
                "end": 1917971,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1917962
            },
            "end": 1917932,
            "length": 311,
            "sequence": {
                "nucleotide": "GTATCCAGTGTGAGGTGAAGCTGGTGGAGTCTGGAGGAGGCTTGGTACAGCCTGGGGGTTCTCTGAGTCTCTCCTGTGCAGCTTCTGGATTCACCTTCACTGATTACTACATGAGCTGGGTCCGCCAGCCTCCAGGGAAGGCACTTGAGTGGTTGGGTTTTATTAGAAACAAAGCTAATGGTTACACAACAGAGTACAGTGCATCTGTGAAGGGTCGGTTCACCATCTCCAGAGATAATTCCCAAAGCATCCTCTATCTTCAAATGAATGCCCTGAGAGCTGAGGACAGTGCCACTTATTACTGTGCAAGATATA",
                "read frame": 1,
                "strand": True
            },
            "start": 1917617
        },
        "accession": "BN000872",
        "allele": "S107.3.62*01",
        "confirmed": True,
        "end": 1917932,
        "family": "S107",
        "framed": True,
        "functionality": "F",
        "gene": "S107.3.62",
        "id": "S107.3.62*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV7-3*01",
        "imgt family name": "IGHV7",
        "imgt gene name": "IGHV7-3",
        "length": 315,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCAGTGTGAGGTGAAGCTGGTGGAGTCTGGAGGAGGCTTGGTACAGCCTGGGGGTTCTCTGAGTCTCTCCTGTGCAGCTTCTGGATTCACCTTCACTGATTACTACATGAGCTGGGTCCGCCAGCCTCCAGGGAAGGCACTTGAGTGGTTGGGTTTTATTAGAAACAAAGCTAATGGTTACACAACAGAGTACAGTGCATCTGTGAAGGGTCGGTTCACCATCTCCAGAGATAATTCCCAAAGCATCCTCTATCTTCAAATGAATGCCCTGAGAGCTGAGGACAGTGCCACTTATTACTGTGCAAGATATA",
            "read frame": 1,
            "strand": True
        },
        "start": 1917617,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1848331,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGCG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1848324
            },
            "3 nonamer": {
                "end": 1848363,
                "length": 9,
                "sequence": {
                    "nucleotide": "AAAAAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1848354
            },
            "end": 1848324,
            "length": 311,
            "sequence": {
                "nucleotide": "GTATCCAGTGTGAGGTGAAGCTGATGGAATCTGGAGGAGGCTTGGTACAGCCTGGGGCTTCTCTGAGACTCTCCTGTGCAGCTTCTGGATTCACCTTTACTGATTACTACATGAGCTGGGTCCGCCAGCCTCCAGGGAAGGCACCTGAGTGGTTGGCTTTGATTAGAAACAAAGCTAATGGTTACACAACAGAGTATACTGCATCTGTTAAGGGTCGGTTCACCATCTCCAGAGATAATTCCCAAAACATCCTCTATCTTCAAATGAACACCCTGAGGGCTGAGGACAGTGCCACTTATTACTGTGTAAAAGCTGTA",
                "read frame": 2,
                "strand": True
            },
            "start": 1848007
        },
        "accession": "BN000872",
        "allele": "S107.4.65*01",
        "confirmed": True,
        "end": 1848324,
        "family": "S107",
        "framed": True,
        "functionality": "F",
        "gene": "S107.4.65",
        "id": "S107.4.65*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV7-4*01",
        "imgt family name": "IGHV7",
        "imgt gene name": "IGHV7-4",
        "length": 317,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCAGTGTGAGGTGAAGCTGATGGAATCTGGAGGAGGCTTGGTACAGCCTGGGGCTTCTCTGAGACTCTCCTGTGCAGCTTCTGGATTCACCTTTACTGATTACTACATGAGCTGGGTCCGCCAGCCTCCAGGGAAGGCACCTGAGTGGTTGGCTTTGATTAGAAACAAAGCTAATGGTTACACAACAGAGTATACTGCATCTGTTAAGGGTCGGTTCACCATCTCCAGAGATAATTCCCAAAACATCCTCTATCTTCAAATGAACACCCTGAGGGCTGAGGACAGTGCCACTTATTACTGTGTAAAAGCTGTA",
            "read frame": 2,
            "strand": True
        },
        "start": 1848007,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2139203,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2139196
            },
            "3 nonamer": {
                "end": 2139234,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCATAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2139225
            },
            "end": 2139196,
            "length": 305,
            "sequence": {
                "nucleotide": "GGGTCAATTCAGAGGTTCAGCTGCAGCAGTCTGGGGCAGAGCTTGTGAGGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTCAACATTAAAGACTACTATATGCACTGGGTGAAGCAGAGGCCTGAACAGGGCCTGGAGTGGATTGGAAGGATTGATCCTGAGGATGGTGATACTGAATATGCCCCGAAGTTCCAGGGCAAGGCCACTATGACTGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTACTACA",
                "read frame": 1,
                "strand": True
            },
            "start": 2138891
        },
        "accession": "BN000872",
        "allele": "SM7.1.44*01",
        "confirmed": True,
        "end": 2139196,
        "family": "SM7",
        "framed": True,
        "functionality": "F",
        "gene": "SM7.1.44",
        "id": "SM7.1.44*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV14-1*01",
        "imgt family name": "IGHV14",
        "imgt gene name": "IGHV14-1",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGGTCAATTCAGAGGTTCAGCTGCAGCAGTCTGGGGCAGAGCTTGTGAGGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTCAACATTAAAGACTACTATATGCACTGGGTGAAGCAGAGGCCTGAACAGGGCCTGGAGTGGATTGGAAGGATTGATCCTGAGGATGGTGATACTGAATATGCCCCGAAGTTCCAGGGCAAGGCCACTATGACTGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTACTACA",
            "read frame": 1,
            "strand": True
        },
        "start": 2138891,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2076650,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2076643
            },
            "3 nonamer": {
                "end": 2076681,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2076672
            },
            "end": 2076643,
            "length": 305,
            "sequence": {
                "nucleotide": "GGGTCAATTCAGAGGTTCAGCTGCAGCAGTCTGGGGCAGAGCTTGTGAAGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTCAACATTAAAGACTACTATATGCACTGGGTGAAGCAGAGGACTGAACAGGGCCTGGAGTGGATTGGAAGGATTGATCCTGAGGATGGTGAAACTAAATATGCCCCGAAATTCCAGGGCAAGGCCACTATAACAGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTGCTAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 2076338
        },
        "accession": "BN000872",
        "allele": "SM7.2.49*01",
        "confirmed": True,
        "end": 2076643,
        "family": "SM7",
        "framed": True,
        "functionality": "F",
        "gene": "SM7.2.49",
        "id": "SM7.2.49*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV14-2*01",
        "imgt family name": "IGHV14",
        "imgt gene name": "IGHV14-2",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGGTCAATTCAGAGGTTCAGCTGCAGCAGTCTGGGGCAGAGCTTGTGAAGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTCAACATTAAAGACTACTATATGCACTGGGTGAAGCAGAGGACTGAACAGGGCCTGGAGTGGATTGGAAGGATTGATCCTGAGGATGGTGAAACTAAATATGCCCCGAAATTCCAGGGCAAGGCCACTATAACAGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTGCTAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 2076338,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2011274,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2011267
            },
            "3 nonamer": {
                "end": 2011305,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2011296
            },
            "end": 2011267,
            "length": 305,
            "sequence": {
                "nucleotide": "GGGTCAATTCAGAGGTTCAGCTGCAGCAGTCTGTGGCAGAGCTTGTGAGGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTCAACATTAAAAACACCTATATGCACTGGGTGAAGCAGAGGCCTGAACAGGGCCTGGAGTGGATTGGAAGGATTGATCCTGCGAATGGTAATACTAAATATGCCCCGAAGTTCCAGGGCAAGGCCACTATAACTGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCATCTATTACTGTGCTAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 2010962
        },
        "accession": "BN000872",
        "allele": "SM7.3.54*01",
        "confirmed": True,
        "end": 2011267,
        "family": "SM7",
        "framed": True,
        "functionality": "F",
        "gene": "SM7.3.54",
        "id": "SM7.3.54*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV14-3*01",
        "imgt family name": "IGHV14",
        "imgt gene name": "IGHV14-3",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGGTCAATTCAGAGGTTCAGCTGCAGCAGTCTGTGGCAGAGCTTGTGAGGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTCAACATTAAAAACACCTATATGCACTGGGTGAAGCAGAGGCCTGAACAGGGCCTGGAGTGGATTGGAAGGATTGATCCTGCGAATGGTAATACTAAATATGCCCCGAAGTTCCAGGGCAAGGCCACTATAACTGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCATCTATTACTGTGCTAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 2010962,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1894681,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1894674
            },
            "3 nonamer": {
                "end": 1894712,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1894703
            },
            "end": 1894674,
            "length": 291,
            "sequence": {
                "nucleotide": "GTTCAGCTGCAGCAGTCTGGGGCTGAGCTTGTGAGGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTTAACATTAAAGACGACTATATGCACTGGGTGAAGCAGAGGCCTGAACAGGGCCTGGAGTGGATTGGATGGATTGATCCTGAGAATGGTGATACTGAATATGCCTCGAAGTTCCAGGGCAAGGCCACTATAACAGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTACTACA",
                "read frame": 0,
                "strand": True
            },
            "start": 1894383
        },
        "accession": "BN000872",
        "allele": "SM7.4.63*01",
        "confirmed": True,
        "end": 1894674,
        "family": "SM7",
        "framed": True,
        "functionality": "F",
        "gene": "SM7.4.63",
        "id": "SM7.4.63*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV14-4*01",
        "imgt family name": "IGHV14",
        "imgt gene name": "IGHV14-4",
        "length": 291,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTTCAGCTGCAGCAGTCTGGGGCTGAGCTTGTGAGGCCAGGGGCCTCAGTCAAGTTGTCCTGCACAGCTTCTGGCTTTAACATTAAAGACGACTATATGCACTGGGTGAAGCAGAGGCCTGAACAGGGCCTGGAGTGGATTGGATGGATTGATCCTGAGAATGGTGATACTGAATATGCCTCGAAGTTCCAGGGCAAGGCCACTATAACAGCAGACACATCCTCCAACACAGCCTACCTGCAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTACTACA",
            "read frame": 0,
            "strand": True
        },
        "start": 1894383,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1977191,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1977184
            },
            "3 nonamer": {
                "end": 1977222,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1977213
            },
            "end": 1977184,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGCCCAAGCACAGATCCAGTTGGTGCAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGGGTATACCTTCACAGAATATCCAATGCACTGGGTGAAGCAGGCTCCAGGAAAGGGTTTCAAGTGGATGGGCATGATATACACCGACACTGGAGAGCCAACATATGCTGAAGAGTTCAAGGGACGGTTTGCCTTCTCTTTGGAGACCTCTGCCAGCACTGCCTATTTGCAGATCAACAACCTCAAAAATGAGGACACGGCTACATATTTCTGTGTAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1976879
        },
        "accession": "BN000872",
        "allele": "VGAM3.8-1-57*01",
        "confirmed": True,
        "end": 1977184,
        "family": "VGAM3.8",
        "framed": True,
        "functionality": "F",
        "gene": "VGAM3.8-1-57",
        "id": "VGAM3.8-1-57*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV9-1*01",
        "imgt family name": "IGHV9",
        "imgt gene name": "IGHV9-1",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGCCCAAGCACAGATCCAGTTGGTGCAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGGGTATACCTTCACAGAATATCCAATGCACTGGGTGAAGCAGGCTCCAGGAAAGGGTTTCAAGTGGATGGGCATGATATACACCGACACTGGAGAGCCAACATATGCTGAAGAGTTCAAGGGACGGTTTGCCTTCTCTTTGGAGACCTCTGCCAGCACTGCCTATTTGCAGATCAACAACCTCAAAAATGAGGACACGGCTACATATTTCTGTGTAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1976879,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1962118,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1962111
            },
            "3 nonamer": {
                "end": 1962149,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1962140
            },
            "end": 1962111,
            "length": 305,
            "sequence": {
                "nucleotide": "GGGCTCAAGCACAGATCCAGTTCGTGCAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGTGTATACCTTCACAGAATATCCAATGCACTGGGTGAAGCAGGCTCCAGGAAAGGGTTTCAAGTGGATGGGCTGGATAAACACCTACTCTGGAGAGCCAACATATGCTGACGACTTCAAGGGACGGTTTGCCTTTTCTTTGGAAACCTCTGCCAGCACTGCCTATTTGCAGATCAACAACCTCAAAAATGAGGACACGGCTACATATTTCTGTGCAAGA",
                "read frame": 2,
                "strand": True
            },
            "start": 1961806
        },
        "accession": "BN000872",
        "allele": "VGAM3.8-2-59*01",
        "confirmed": True,
        "end": 1962111,
        "family": "VGAM3.8",
        "framed": True,
        "functionality": "F",
        "gene": "VGAM3.8-2-59",
        "id": "VGAM3.8-2-59*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV9-2*01",
        "imgt family name": "IGHV9",
        "imgt gene name": "IGHV9-2",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGGCTCAAGCACAGATCCAGTTCGTGCAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGTGTATACCTTCACAGAATATCCAATGCACTGGGTGAAGCAGGCTCCAGGAAAGGGTTTCAAGTGGATGGGCTGGATAAACACCTACTCTGGAGAGCCAACATATGCTGACGACTTCAAGGGACGGTTTGCCTTTTCTTTGGAAACCTCTGCCAGCACTGCCTATTTGCAGATCAACAACCTCAAAAATGAGGACACGGCTACATATTTCTGTGCAAGA",
            "read frame": 2,
            "strand": True
        },
        "start": 1961806,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1930427,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1930420
            },
            "3 nonamer": {
                "end": 1930458,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAAAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1930449
            },
            "end": 1930420,
            "length": 305,
            "sequence": {
                "nucleotide": "GTGCCCAAGCACAGATCCAGTTGGTACAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGGGTATACCTTCACAACCTATGGAATGAGCTGGGTGAAACAGGCTCCAGGAAAGGGTTTAAAGTGGATGGGCTGGATAAACACCTACTCTGGAGTGCCAACATATGCTGATGACTTCAAGGGACGGTTTGCCTTCTCTTTGGAAACCTCTGCCAGCACTGCCTATTTGCAGATCAACAACCTCAAAAATGAGGACACGGCTACATATTTCTGTGCAAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1930115
        },
        "accession": "BN000872",
        "allele": "VGAM3.8-3-61*01",
        "confirmed": True,
        "end": 1930420,
        "family": "VGAM3.8",
        "framed": True,
        "functionality": "F",
        "gene": "VGAM3.8-3-61",
        "id": "VGAM3.8-3-61*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV9-3*01",
        "imgt family name": "IGHV9",
        "imgt gene name": "IGHV9-3",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGCCCAAGCACAGATCCAGTTGGTACAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGGGTATACCTTCACAACCTATGGAATGAGCTGGGTGAAACAGGCTCCAGGAAAGGGTTTAAAGTGGATGGGCTGGATAAACACCTACTCTGGAGTGCCAACATATGCTGATGACTTCAAGGGACGGTTTGCCTTCTCTTTGGAAACCTCTGCCAGCACTGCCTATTTGCAGATCAACAACCTCAAAAATGAGGACACGGCTACATATTTCTGTGCAAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1930115,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1771155,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1771148
            },
            "3 nonamer": {
                "end": 1771186,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAGAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1771177
            },
            "end": 1771148,
            "length": 305,
            "sequence": {
                "nucleotide": "GTATCCAAGCACAGATCCAGTTGGTGCAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGGGTATACCTTCACAACTGCTGGAATGCAGTGGGTGCAAAAGATGCCAGGAAAGGGTTTTAAGTGGATTGGCTGGATAAACACCCACTCTGGAGAGCCAAAATATGCAGAAGACTTCAAGGGACGGTTTGCCTTCTCTTTGGAAACCTCTGCCAGCACTGCCTATTTACAGATAAGCAACCTCAAAAATGAGGACACGGCTACGTATTTCTGTGCGAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1770843
        },
        "accession": "BN000872",
        "allele": "VGAM3.8-4-71*01",
        "confirmed": True,
        "end": 1771148,
        "family": "VGAM3.8",
        "framed": True,
        "functionality": "F",
        "gene": "VGAM3.8-4-71",
        "id": "VGAM3.8-4-71*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV9-4*01",
        "imgt family name": "IGHV9",
        "imgt gene name": "IGHV9-4",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCAAGCACAGATCCAGTTGGTGCAGTCTGGACCTGAGCTGAAGAAGCCTGGAGAGACAGTCAAGATCTCCTGCAAGGCTTCTGGGTATACCTTCACAACTGCTGGAATGCAGTGGGTGCAAAAGATGCCAGGAAAGGGTTTTAAGTGGATTGGCTGGATAAACACCCACTCTGGAGAGCCAAAATATGCAGAAGACTTCAAGGGACGGTTTGCCTTCTCTTTGGAAACCTCTGCCAGCACTGCCTATTTACAGATAAGCAACCTCAAAAATGAGGACACGGCTACGTATTTCTGTGCGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1770843,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1592111,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1592104
            },
            "3 nonamer": {
                "end": 1592143,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1592134
            },
            "end": 1592104,
            "length": 313,
            "sequence": {
                "nucleotide": "GTGTGCATTGTGAGGTGCAGCTTGTTGAGTCTGGTGGAGGATTGGTGCAGCCTAAAGGGTCATTGAAACTCTCATGTGCAGCCTCTGGATTCAGCTTCAATACCTACGCCATGAACTGGGTCCGCCAGGCTCCAGGAAAGGGTTTGGAATGGGTTGCTCGCATAAGAAGTAAAAGTAATAATTATGCAACATATTATGCCGATTCAGTGAAAGACAGATTCACCATCTCCAGAGATGATTCAGAAAGCATGCTCTATCTGCAAATGAACAACTTGAAAACTGAGGACACAGCCATGTATTACTGTGTGAGACA",
                "read frame": 0,
                "strand": True
            },
            "start": 1591791
        },
        "accession": "BN000872",
        "allele": "VH10.1.86*01",
        "confirmed": True,
        "end": 1592104,
        "family": "VH10",
        "framed": True,
        "functionality": "F",
        "gene": "VH10.1.86",
        "id": "VH10.1.86*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV10-1*01",
        "imgt family name": "IGHV10",
        "imgt gene name": "IGHV10-1",
        "length": 313,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTGCATTGTGAGGTGCAGCTTGTTGAGTCTGGTGGAGGATTGGTGCAGCCTAAAGGGTCATTGAAACTCTCATGTGCAGCCTCTGGATTCAGCTTCAATACCTACGCCATGAACTGGGTCCGCCAGGCTCCAGGAAAGGGTTTGGAATGGGTTGCTCGCATAAGAAGTAAAAGTAATAATTATGCAACATATTATGCCGATTCAGTGAAAGACAGATTCACCATCTCCAGAGATGATTCAGAAAGCATGCTCTATCTGCAAATGAACAACTTGAAAACTGAGGACACAGCCATGTATTACTGTGTGAGACA",
            "read frame": 0,
            "strand": True
        },
        "start": 1591791,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1575113,
                "length": 7,
                "sequence": {
                    "nucleotide": "TATAGGG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1575106
            },
            "3 nonamer": {
                "end": 1575145,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1575136
            },
            "end": 1575106,
            "length": 301,
            "sequence": {
                "nucleotide": "GTGTGCATTGTGAGGTACAGCTTGTTGGGTCTGGTGAAGGATTGGTGCAGCCTAAAGGGTCACTTGTATTCAGCCTCTGGATTCATCTTCAATACATACACCTTGGAGTGGGTGTGTCAGACTCCAAGAAAGGGTCTGGAATGGGTTGCATGCATAAGAACAAAAAGTAATAATTATGCCACATATACTGATTCAGTGAAAGACAGATTCACCATCTCCAGAGATGATTCTCAAAGCATGGTCAACCTGCAAATGAACAATTTGAAAACTGAGGACATAGACCTGTATTAATTACAAGAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1574805
        },
        "accession": "BN000872",
        "allele": "VH10.2pg.89*01",
        "confirmed": True,
        "end": 1575106,
        "family": "VH10",
        "framed": False,
        "functionality": "P",
        "gene": "VH10.2pg.89",
        "id": "VH10.2pg.89*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV10-2*01",
        "imgt family name": "IGHV10",
        "imgt gene name": "IGHV10-2",
        "length": 301,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTGCATTGTGAGGTACAGCTTGTTGGGTCTGGTGAAGGATTGGTGCAGCCTAAAGGGTCACTTGTATTCAGCCTCTGGATTCATCTTCAATACATACACCTTGGAGTGGGTGTGTCAGACTCCAAGAAAGGGTCTGGAATGGGTTGCATGCATAAGAACAAAAAGTAATAATTATGCCACATATACTGATTCAGTGAAAGACAGATTCACCATCTCCAGAGATGATTCTCAAAGCATGGTCAACCTGCAAATGAACAATTTGAAAACTGAGGACATAGACCTGTATTAATTACAAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1574805,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1547675,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1547668
            },
            "3 nonamer": {
                "end": 1547707,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1547698
            },
            "end": 1547668,
            "length": 313,
            "sequence": {
                "nucleotide": "GTGTGCATTGTGAGGTGCAGCTTGTTGAGTCTGGTGGAGGATTGGTGCAGCCTAAAGGATCATTGAAACTCTCATGTGCCGCCTCTGGTTTCACCTTCAATACCTATGCCATGCACTGGGTCCGCCAGGCTCCAGGAAAGGGTTTGGAATGGGTTGCTCGCATAAGAAGTAAAAGTAGTAATTATGCAACATATTATGCCGATTCAGTGAAAGACAGATTCACCATCTCCAGAGATGATTCACAAAGCATGCTCTATCTGCAAATGAACAACCTGAAAACTGAGGACACAGCCATGTATTACTGTGTGAGAGA",
                "read frame": 0,
                "strand": True
            },
            "start": 1547355
        },
        "accession": "BN000872",
        "allele": "VH10.3.91*01",
        "confirmed": True,
        "end": 1547668,
        "family": "VH10",
        "framed": True,
        "functionality": "F",
        "gene": "VH10.3.91",
        "id": "VH10.3.91*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV10-3*01",
        "imgt family name": "IGHV10",
        "imgt gene name": "IGHV10-3",
        "length": 313,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTGCATTGTGAGGTGCAGCTTGTTGAGTCTGGTGGAGGATTGGTGCAGCCTAAAGGATCATTGAAACTCTCATGTGCCGCCTCTGGTTTCACCTTCAATACCTATGCCATGCACTGGGTCCGCCAGGCTCCAGGAAAGGGTTTGGAATGGGTTGCTCGCATAAGAAGTAAAAGTAGTAATTATGCAACATATTATGCCGATTCAGTGAAAGACAGATTCACCATCTCCAGAGATGATTCACAAAGCATGCTCTATCTGCAAATGAACAACCTGAAAACTGAGGACACAGCCATGTATTACTGTGTGAGAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1547355,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2089240,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAATG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2089233
            },
            "3 nonamer": {
                "end": 2089272,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACATAAACC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2089263
            },
            "end": 2089233,
            "length": 307,
            "sequence": {
                "nucleotide": "ATGTCCAGTGTGAAGTGCAGCTGTTGGAGACTGGAGAAGGCTTGGTGCCACCTGGGGGGTCACGGGGACTCTCTTGTGAAGGCTCAGGGTTCACTTTTAGTGGCTTCTGGATGAGCTGGGTTCGACAGACACCTGGGAAGACCCTGGAGTGGATTGGAGACATTAATTCTGATGGCAGTGCAATAAACTACGCACCATCCATAAAGGATCGCTTCACTATCTTCAGAGACAATGACAAGAGCACCCTGTACCTGCAGATGAGCAATGTGCGATCTGAGGACACAGCCACGTATTTCTGTATGAGATA",
                "read frame": 1,
                "strand": True
            },
            "start": 2088926
        },
        "accession": "BN000872",
        "allele": "VH11.1.48*01",
        "confirmed": True,
        "end": 2089233,
        "family": "VH11",
        "framed": True,
        "functionality": "F",
        "gene": "VH11.1.48",
        "id": "VH11.1.48*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV11-1*01",
        "imgt family name": "IGHV11",
        "imgt gene name": "IGHV11-1",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATGTCCAGTGTGAAGTGCAGCTGTTGGAGACTGGAGAAGGCTTGGTGCCACCTGGGGGGTCACGGGGACTCTCTTGTGAAGGCTCAGGGTTCACTTTTAGTGGCTTCTGGATGAGCTGGGTTCGACAGACACCTGGGAAGACCCTGGAGTGGATTGGAGACATTAATTCTGATGGCAGTGCAATAAACTACGCACCATCCATAAAGGATCGCTTCACTATCTTCAGAGACAATGACAAGAGCACCCTGTACCTGCAGATGAGCAATGTGCGATCTGAGGACACAGCCACGTATTTCTGTATGAGATA",
            "read frame": 1,
            "strand": True
        },
        "start": 2088926,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2022878,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAATG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2022871
            },
            "3 nonamer": {
                "end": 2022910,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACATAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2022901
            },
            "end": 2022871,
            "length": 307,
            "sequence": {
                "nucleotide": "ATGTCCAGTGTGAAGTGCAGCTGTTGGAGACTGGAGGAGGCTTGGTGCAACCTGGGGGGTCACGGGGACTCTCTTGTGAAGGCTCAGGGTTCACTTTTAGTGGCTTCTGGATGAGCTGGGTTCGACAGACACCTGGGAAGACCCTGGAGTGGATTGGAGACATTAATTCTGATGGCAGTGCAATAAACTACGCACCATCCATAAAGGATCGATTCACTATCTTCAGAGACAATGACAAGAGCACCCTGTACCTGCAGATGAGCAATGTGCGATCGGAGGACACAGCCACGTATTTCTGTATGAGATA",
                "read frame": 0,
                "strand": True
            },
            "start": 2022564
        },
        "accession": "BN000872",
        "allele": "VH11.2.53*01",
        "confirmed": True,
        "end": 2022871,
        "family": "VH11",
        "framed": True,
        "functionality": "F",
        "gene": "VH11.2.53",
        "id": "VH11.2.53*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV11-2*01",
        "imgt family name": "IGHV11",
        "imgt gene name": "IGHV11-2",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "ATGTCCAGTGTGAAGTGCAGCTGTTGGAGACTGGAGGAGGCTTGGTGCAACCTGGGGGGTCACGGGGACTCTCTTGTGAAGGCTCAGGGTTCACTTTTAGTGGCTTCTGGATGAGCTGGGTTCGACAGACACCTGGGAAGACCCTGGAGTGGATTGGAGACATTAATTCTGATGGCAGTGCAATAAACTACGCACCATCCATAAAGGATCGATTCACTATCTTCAGAGACAATGACAAGAGCACCCTGTACCTGCAGATGAGCAATGTGCGATCGGAGGACACAGCCACGTATTTCTGTATGAGATA",
            "read frame": 0,
            "strand": True
        },
        "start": 2022564,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1963860,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1963853
            },
            "3 nonamer": {
                "end": 1963892,
                "length": 9,
                "sequence": {
                    "nucleotide": "CATAAACTT",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1963883
            },
            "end": 1963853,
            "length": 316,
            "sequence": {
                "nucleotide": "GTGTTCTGTCCAAGATTTAGCTTAAGGAGTCAGGATCTGCTCTCATCAAGCCATCACAGCCACTGTCTCTCACCTGCACAGTCTCTGGATTCTCCATCACAAGTAGTAGTTATTGCTGGCACTGGATCTGCCAGCCCCCAGGAAAGGGGTTAGAATGGATGGGGCACATATGTTATGAAGGTTCACTAAACTATAGTCCATCCCTCAAAAGCCGCAGCACCATCTCCAGAGACACATCTCTGAACAAATTCTTTATCCAGCTGAGCTCTCTGACTGATGAGGACACAGTCATGTACTACTGTTCTAGGGAAAACCA",
                "read frame": 2,
                "strand": True
            },
            "start": 1963537
        },
        "accession": "BN000872",
        "allele": "PG.11.58*01",
        "confirmed": True,
        "end": 1963853,
        "family": "VH12",
        "framed": True,
        "functionality": "P",
        "gene": "PG.11.58",
        "id": "PG.11.58*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV12-1*01",
        "imgt family name": "IGHV12",
        "imgt gene name": "IGHV12-1",
        "length": 316,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTGTTCTGTCCAAGATTTAGCTTAAGGAGTCAGGATCTGCTCTCATCAAGCCATCACAGCCACTGTCTCTCACCTGCACAGTCTCTGGATTCTCCATCACAAGTAGTAGTTATTGCTGGCACTGGATCTGCCAGCCCCCAGGAAAGGGGTTAGAATGGATGGGGCACATATGTTATGAAGGTTCACTAAACTATAGTCCATCCCTCAAAAGCCGCAGCACCATCTCCAGAGACACATCTCTGAACAAATTCTTTATCCAGCTGAGCTCTCTGACTGATGAGGACACAGTCATGTACTACTGTTCTAGGGAAAACCA",
            "read frame": 2,
            "strand": True
        },
        "start": 1963537,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 1943571,
            "length": 306,
            "sequence": {
                "nucleotide": "ATTCTGTCCCAGATTCAGCTTAAGGAATCAGGACCTGCTGTCATCAAGCCATCATAGTCACTGTCTCTCACCTGCACAGTCTCTGGATTCTCCATCACTAGTAGTGGTTTTTGCTGGCACTGGATATGCCAGCCCCCAGGAAAGGGGTTAGAGTGGATGGGGCGCATATGTTATGAAGGTTCCATATACTATAGTCCATCCCTCAAAAGTCCCAGCACCATCTCCAGAGACATATCACTGAATAAATTCTTTATCCAGCTGAGCTCTGTAACTGATGAGGACACAGCCATGTACTACTATTCCAGG",
                "read frame": 0,
                "strand": True
            },
            "start": 1943265
        },
        "accession": "BN000872",
        "allele": "PG.12.60*01",
        "confirmed": True,
        "end": 1943579,
        "family": "VH12",
        "framed": True,
        "functionality": "P",
        "gene": "PG.12.60",
        "id": "PG.12.60*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV12-2*01",
        "imgt family name": "IGHV12",
        "imgt gene name": "IGHV12-2",
        "length": 305,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CAGATTCAGCTTAAGGAATCAGGACCTGCTGTCATCAAGCCATCATAGTCACTGTCTCTCACCTGCACAGTCTCTGGATTCTCCATCACTAGTAGTGGTTTTTGCTGGCACTGGATATGCCAGCCCCCAGGAAAGGGGTTAGAGTGGATGGGGCGCATATGTTATGAAGGTTCCATATACTATAGTCCATCCCTCAAAAGTCCCAGCACCATCTCCAGAGACATATCACTGAATAAATTCTTTATCCAGCTGAGCTCTGTAACTGATGAGGACACAGCCATGTACTACTATTCCAGGGAAAACCA",
            "read frame": 0,
            "strand": True
        },
        "start": 1943274,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1704596,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAATG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1704589
            },
            "3 nonamer": {
                "end": 1704628,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACACAAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1704619
            },
            "end": 1704589,
            "length": 311,
            "sequence": {
                "nucleotide": "GTAGCCTGTCTCAGATGCAGCTTCAGGAGTCAGGACCTGGCCTGGTGAAACCCTCACAGTCACTCTTCCTTACCTGCTCTATTACTGGTTTCCCCATCACCAGTGGTTACTACTGGATCTGGATCCGTCAGTCACCTGGGAAACCCCTAGAATGGATGGGGTACATCACTCATAGTGGGGAAACTTTCTACAACCCATCTCTCCAGAGCCCCATCTCCATTACTAGAGAAACGTCAAAGAACCAGTTCTTCCTCCAATTGAACTCTGTGACCACAGAGGACACAGCCATGTATTACTGTGCAGGAGACAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 1704278
        },
        "accession": "BN000872",
        "allele": "VH12.1.78*01",
        "confirmed": True,
        "end": 1704589,
        "family": "VH12",
        "framed": True,
        "functionality": "F",
        "gene": "VH12.1.78",
        "id": "VH12.1.78*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV12-3*01",
        "imgt family name": "IGHV12",
        "imgt gene name": "IGHV12-3",
        "length": 311,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTAGCCTGTCTCAGATGCAGCTTCAGGAGTCAGGACCTGGCCTGGTGAAACCCTCACAGTCACTCTTCCTTACCTGCTCTATTACTGGTTTCCCCATCACCAGTGGTTACTACTGGATCTGGATCCGTCAGTCACCTGGGAAACCCCTAGAATGGATGGGGTACATCACTCATAGTGGGGAAACTTTCTACAACCCATCTCTCCAGAGCCCCATCTCCATTACTAGAGAAACGTCAAAGAACCAGTTCTTCCTCCAATTGAACTCTGTGACCACAGAGGACACAGCCATGTATTACTGTGCAGGAGACAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 1704278,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "accession": "BN000872",
        "allele": "IGHV15-1*01",
        "confirmed": True,
        "end": 1913087,
        "family": "VH15",
        "framed": False,
        "functionality": "P",
        "gene": "IGHV15-1",
        "id": "IGHV15-1*01*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV15-1*01",
        "imgt family name": "IGHV15",
        "imgt gene name": "IGHV15-1",
        "length": 291,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "CACTCCATTGTTTAGCTGCAGCAGTCTGAAGCTGCGCTGAGGACTCTTGGAGCTTCAGAATAGGTGCCCTTCAAATGCTGTGATATGGGTAGCTTTCCCTTTTCCTTTAGGAATTGCATGACATAGAACCCTGAATGTAATGGAAACATAGACCCAAGCATGAGAAATACACTCTATGGACAGAATTAACAAGGCAGAGTCACAATGGATGCAGACAAAATGTCCAAAACTGCCTACATGTAGGTCAATAATCTGACAGCTGAGGACTCTTCTATCACTGCAGAAGGAAGA",
            "read frame": 0,
            "strand": True
        },
        "start": 1912796,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 1506537,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1506530
            },
            "3 nonamer": {
                "end": 1506569,
                "length": 9,
                "sequence": {
                    "nucleotide": "TCAAAAAGT",
                    "read frame": 1,
                    "strand": True
                },
                "source": "BN000872",
                "start": 1506560
            },
            "end": 1506530,
            "length": 308,
            "sequence": {
                "nucleotide": "GTATCCAATCCCAGGTTCACCTACAACAGTCTGGTTCTGAACTGAGGAGTCCTGGGTCTTCAGTAAAGCTGTCATGCAAGGATTTTGATTCAGAAGTCTTCCCTATTGCTTATATGAGTTGGGTAAGGCAGAAGCCTGGGCATGGATTTGAATGGATTGGAGGCATACTCCCAAGTATTGGTAGAACAATCTATGGAGAGAAGTTTGAGGACAAAGCCACATTGGATGCAGACACACTGTCCAACACAGCCTACTTGGAGCTCAACAGTCTGACATCCGAGGACTCTGCTATCTACTACTGTGCAAGG",
                "read frame": 0,
                "strand": True
            },
            "start": 1506222
        },
        "accession": "BN000872",
        "allele": "VH15.1.95*01",
        "confirmed": True,
        "end": 1506530,
        "family": "VH15",
        "framed": True,
        "functionality": "F",
        "gene": "VH15.1.95",
        "id": "VH15.1.95*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV15-2*01",
        "imgt family name": "IGHV15",
        "imgt gene name": "IGHV15-2",
        "length": 308,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GTATCCAATCCCAGGTTCACCTACAACAGTCTGGTTCTGAACTGAGGAGTCCTGGGTCTTCAGTAAAGCTGTCATGCAAGGATTTTGATTCAGAAGTCTTCCCTATTGCTTATATGAGTTGGGTAAGGCAGAAGCCTGGGCATGGATTTGAATGGATTGGAGGCATACTCCCAAGTATTGGTAGAACAATCTATGGAGAGAAGTTTGAGGACAAAGCCACATTGGATGCAGACACACTGTCCAACACAGCCTACTTGGAGCTCAACAGTCTGACATCCGAGGACTCTGCTATCTACTACTGTGCAAGG",
            "read frame": 0,
            "strand": True
        },
        "start": 1506222,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2002291,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2002284
            },
            "3 nonamer": {
                "end": 2002323,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACAAAACCC",
                    "read frame": 0,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2002314
            },
            "end": 2002284,
            "length": 310,
            "sequence": {
                "nucleotide": "GCATACAGGGTGAGGTGCAGCTGGTGGAATCTGGAGGCAGCTTGGGACAGCCTGGAGGGTCCACTAAACTCTCTTGTGAAGAAGCATCTGGATTCACTTTCAGTGATCATTGGATGGACTGGTTTCGCCAAGCCCCAGGCATGAGGCTAGAATGGTTAGCAAATACAAACCATGATGAGAGTGGAAAAGGCTATGCAGAGTCTGTGAAAGACAGATTCTCCATCTCCAGAGACAATTCTGAGAACTTATTGTATCTACAAATGAACAGTCTGAGAAACGAAGACACGGCTCTGTATTATTGTGCCAGAGA",
                "read frame": 1,
                "strand": True
            },
            "start": 2001974
        },
        "accession": "BN000872",
        "allele": "VH16.1.55*01",
        "confirmed": True,
        "end": 2002284,
        "family": "VH16",
        "framed": True,
        "functionality": "F",
        "gene": "VH16.1.55",
        "id": "VH16.1.55*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV16-1*01",
        "imgt family name": "IGHV16",
        "imgt gene name": "IGHV16-1",
        "length": 310,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GCATACAGGGTGAGGTGCAGCTGGTGGAATCTGGAGGCAGCTTGGGACAGCCTGGAGGGTCCACTAAACTCTCTTGTGAAGAAGCATCTGGATTCACTTTCAGTGATCATTGGATGGACTGGTTTCGCCAAGCCCCAGGCATGAGGCTAGAATGGTTAGCAAATACAAACCATGATGAGAGTGGAAAAGGCTATGCAGAGTCTGTGAAAGACAGATTCTCCATCTCCAGAGACAATTCTGAGAACTTATTGTATCTACAAATGAACAGTCTGAGAAACGAAGACACGGCTCTGTATTATTGTGCCAGAGA",
            "read frame": 1,
            "strand": True
        },
        "start": 2001974,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2122874,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2122867
            },
            "3 nonamer": {
                "end": 2122906,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACATGAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2122897
            },
            "end": 2122867,
            "length": 305,
            "sequence": {
                "nucleotide": "GGGTCCAGTGTGAGGTGAAGCTTCTCCAGTCTGGAGGTGGCCTGGTGCAGCCTGGAGGATCCCTGAAACTCTCCTGTGCAGCCTCAGGAATCGATTTTAGTAGATACTGGATGAGTTGGGTTCGGCGGGCTCCAGGGAAAGGACTAGAATGGATTGGAGAAATTAATCCAGATAGCAGTACAATAAACTATGCACCATCTCTAAAGGATAAATTCATCATCTCCAGAGACAACGCCAAAAATACGCTGTACCTGCAAATGAGCAAAGTGAGATCTGAGGACACAGCCCTTTATTACTGTGCAAGACC",
                "read frame": 0,
                "strand": True
            },
            "start": 2122560
        },
        "accession": "BN000872",
        "allele": "X24.1pg.45*01",
        "confirmed": True,
        "end": 2122867,
        "family": "X24",
        "framed": True,
        "functionality": "F",
        "gene": "X24.1pg.45",
        "id": "X24.1pg.45*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV4-1*01",
        "imgt family name": "IGHV4",
        "imgt gene name": "IGHV4-1",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGGTCCAGTGTGAGGTGAAGCTTCTCCAGTCTGGAGGTGGCCTGGTGCAGCCTGGAGGATCCCTGAAACTCTCCTGTGCAGCCTCAGGAATCGATTTTAGTAGATACTGGATGAGTTGGGTTCGGCGGGCTCCAGGGAAAGGACTAGAATGGATTGGAGAAATTAATCCAGATAGCAGTACAATAAACTATGCACCATCTCTAAAGGATAAATTCATCATCTCCAGAGACAACGCCAAAAATACGCTGTACCTGCAAATGAGCAAAGTGAGATCTGAGGACACAGCCCTTTATTACTGTGCAAGACC",
            "read frame": 0,
            "strand": True
        },
        "start": 2122560,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "3 heptamer": {
                "end": 2057975,
                "length": 7,
                "sequence": {
                    "nucleotide": "CACAGTG",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2057968
            },
            "3 nonamer": {
                "end": 2058007,
                "length": 9,
                "sequence": {
                    "nucleotide": "ACATGAACC",
                    "read frame": 2,
                    "strand": True
                },
                "source": "BN000872",
                "start": 2057998
            },
            "end": 2057968,
            "length": 305,
            "sequence": {
                "nucleotide": "GGGTCCAGTGTGAGGTGAAGCTTCTCGAGTCTGGAGGTGGCCTGGTGCAGCCTGGAGGATCCCTGAAACTCTCCTGTGCAGCCTCAGGATTCGATTTTAGTAAAGACTGGATGAGTTGGGTCCGGCAGGCTACAGGGAAAGGGCTAGAATGAATTGGAGAAATTAATCCAGGTAGCAGTACGATAAACTATACTCCATCTCTAAAGGATAAATTCATCATCTCCAGAGACAACGCCAAAAATACGCTGTACCTGCAAATGAGCAAAGTGAGATCTGAGGACACAGCCCTTTATTACTGTGCAAGACC",
                "read frame": 0,
                "strand": True
            },
            "start": 2057661
        },
        "accession": "BN000872",
        "allele": "X24.2.50*01",
        "confirmed": True,
        "end": 2057968,
        "family": "X24",
        "framed": True,
        "functionality": "P",
        "gene": "X24.2.50",
        "id": "X24.2.50*01-C57BL/6",
        "identified": True,
        "imgt allele name": "IGHV4-2*01",
        "imgt family name": "IGHV4",
        "imgt gene name": "IGHV4-2",
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGGTCCAGTGTGAGGTGAAGCTTCTCGAGTCTGGAGGTGGCCTGGTGCAGCCTGGAGGATCCCTGAAACTCTCCTGTGCAGCCTCAGGATTCGATTTTAGTAAAGACTGGATGAGTTGGGTCCGGCAGGCTACAGGGAAAGGGCTAGAATGAATTGGAGAAATTAATCCAGGTAGCAGTACGATAAACTATACTCCATCTCTAAAGGATAAATTCATCATCTCCAGAGACAACGCCAAAAATACGCTGTACCTGCAAATGAGCAAAGTGAGATCTGAGGACACAGCCCTTTATTACTGTGCAAGACC",
            "read frame": 0,
            "strand": True
        },
        "start": 2057661,
        "strain": "C57BL/6",
        "strand": True,
        "verified": True
    },
    {
        "BN000872": {
            "end": 1786429,
            "length": 311,
            "sequence": {
                "nucleotide": "GTGAGGTGCAGGTGGTAGAATCCGGAGGCAGCTTGATTCAGCTGGGGGGTGTGTGTCGATTAATCTCTCTTGTGAAGCTTCTGGATTCACCTTCAGTAATTACTGGATGGACTGGATTTGCCAAGCCTCAAGGAAGGGGTTAGAGTTGTTAGCAAAATCAAGAGGAATTCTATTGTCACAGCATTCTGGTGTTGCCCTGCACAAAAAGGAGAGAGATAAAGACCCCACGAAGCCTGATTGTTGAGTTTTTATACAGTTTTCAGAGCAGCACCCATTAGGAACAATGTTATTGGTAGAACGGTGTGACTTTT",
                "read frame": 1,
                "strand": True
            },
            "start": 1786118
        },
        "accession": "BN000872",
        "allele": "PG.13.69*01",
        "confirmed": True,
        "end": 1786429,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.13.69",
        "id": "PG.13.69*01-C57BL/6",
        "identified": True,
        "length": 310,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGAGGTGCAGGTGGTAGAATCCGGAGGCAGCTTGATTCAGCTGGGGGGTGTGTGTCGATTAATCTCTCTTGTGAAGCTTCTGGATTCACCTTCAGTAATTACTGGATGGACTGGATTTGCCAAGCCTCAAGGAAGGGGTTAGAGTTGTTAGCAAAATCAAGAGGAATTCTATTGTCACAGCATTCTGGTGTTGCCCTGCACAAAAAGGAGAGAGATAAAGACCCCACGAAGCCTGATTGTTGAGTTTTTATACAGTTTTCAGAGCAGCACCCATTAGGAACAATGTTATTGGTAGAACGGTGTGACTTTT",
            "read frame": 0,
            "strand": True
        },
        "start": 1786119,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "end": 658151,
            "length": 196,
            "sequence": {
                "nucleotide": "ATGGGATGGAGCTATATCATCTTCTTCCTTGTAGCAACAGCTATATGTAAGGGTCTCACAGTAGCAGATTTGAGAGGTCTGGCCATACACTCACATGACAATGACATCCATTCTATCCTTCCATCCACAGGTATCCACTCCCAGGTATAGCTGCATCAATCTGGGGCTGAGCTGGTGAAGCCTGGGATCTCTGTTA",
                "read frame": 2,
                "strand": True
            },
            "start": 657955
        },
        "accession": "BN000872",
        "allele": "PG.19.161*01",
        "confirmed": True,
        "end": 658151,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.19.161",
        "id": "PG.19.161*01-C57BL/6",
        "identified": True,
        "length": 195,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TGGGATGGAGCTATATCATCTTCTTCCTTGTAGCAACAGCTATATGTAAGGGTCTCACAGTAGCAGATTTGAGAGGTCTGGCCATACACTCACATGACAATGACATCCATTCTATCCTTCCATCCACAGGTATCCACTCCCAGGTATAGCTGCATCAATCTGGGGCTGAGCTGGTGAAGCCTGGGATCTCTGTTA",
            "read frame": 0,
            "strand": True
        },
        "start": 657956,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "end": 2352251,
            "length": 315,
            "sequence": {
                "nucleotide": "TGGGCTCAGCTGTGTTTTCCTTGTCCTCATTTTAAGAGGTAATTTGTAGAAATAAGATCCTGCCTGTTTTGTGTACAGGAGAAATAGAAATTTTTTTCTTTCCTCTACTGTGTTTTGTTTTGTTAGTGACAGTTTACAAATAAGCATTCTCTGTTGTGAGGTGTCCAGTGTGAGGTGAAGATGGTGGAGTCTGTGGGAGGCTTAGTGCAGCCTGTAGGGTCATTGAAACCCTCCTGTGCAGCCTCGGGATTCATTCTCACTGACTACTGAATGACCTGGATCCTTCAGGCTTCAAAGAAAAGGATGGAGAGGGTG",
                "read frame": 1,
                "strand": True
            },
            "start": 2351936
        },
        "accession": "BN000872",
        "allele": "PG.3.23*01",
        "confirmed": True,
        "end": 2352251,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.3.23",
        "id": "PG.3.23*01-C57BL/6",
        "identified": True,
        "length": 314,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GGGCTCAGCTGTGTTTTCCTTGTCCTCATTTTAAGAGGTAATTTGTAGAAATAAGATCCTGCCTGTTTTGTGTACAGGAGAAATAGAAATTTTTTTCTTTCCTCTACTGTGTTTTGTTTTGTTAGTGACAGTTTACAAATAAGCATTCTCTGTTGTGAGGTGTCCAGTGTGAGGTGAAGATGGTGGAGTCTGTGGGAGGCTTAGTGCAGCCTGTAGGGTCATTGAAACCCTCCTGTGCAGCCTCGGGATTCATTCTCACTGACTACTGAATGACCTGGATCCTTCAGGCTTCAAAGAAAAGGATGGAGAGGGTG",
            "read frame": 0,
            "strand": True
        },
        "start": 2351937,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "end": 2303285,
            "length": 317,
            "sequence": {
                "nucleotide": "AAAAGGTAATTCATAAAGATAAGATTCTGTCTGTTGTGTGCACATGAGAAACAGAAAAATTGTATTGTTTCTCTATTTGGTTTTGTTTTGTTAGTGACAGTTTCTGACTCAGAATTCTCTGTTTGAAGGTGCCCAGTGTGAGGTGAAGCTTGAGAAGTCTAGGGGAGGCTTAGTGCAGCCTGGAAGGTCCGTGATACGCTCATGTGCAGCCTCTTGATGCACTGCCAGTGACGACTGGTTTGCCTGGGTTTGCCAGGCTCCAAAGAAGGGGCTGCAGTGGGGGTGGGGCATAATTTTTCATGGCTGTGGTAGCCCCT",
                "read frame": 0,
                "strand": True
            },
            "start": 2302968
        },
        "accession": "BN000872",
        "allele": "PG.4.28*01",
        "confirmed": True,
        "end": 2303285,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.4.28",
        "id": "PG.4.28*01-C57BL/6",
        "identified": True,
        "length": 316,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "AAAGGTAATTCATAAAGATAAGATTCTGTCTGTTGTGTGCACATGAGAAACAGAAAAATTGTATTGTTTCTCTATTTGGTTTTGTTTTGTTAGTGACAGTTTCTGACTCAGAATTCTCTGTTTGAAGGTGCCCAGTGTGAGGTGAAGCTTGAGAAGTCTAGGGGAGGCTTAGTGCAGCCTGGAAGGTCCGTGATACGCTCATGTGCAGCCTCTTGATGCACTGCCAGTGACGACTGGTTTGCCTGGGTTTGCCAGGCTCCAAAGAAGGGGCTGCAGTGGGGGTGGGGCATAATTTTTCATGGCTGTGGTAGCCCCT",
            "read frame": 0,
            "strand": True
        },
        "start": 2302969,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "end": 2299106,
            "length": 305,
            "sequence": {
                "nucleotide": "AGACCACTCACCATAGACTTTGGGCTCAGCTTTGCTTTCCTTGTTCTTATTTTAAAAGGTAATTCTTAGAAATAAGATCCTGCCTGTTTTGTGTACATGAGAAATAGAAAAATTTGTTTTCTTTCCTCTATTTTGTTTTGTTTTGTTTTGTTAGTGACAGTTTCCAAATCAGCATTCTCTGCTTTGAGGTGCCCAGTGTGAGGTGAAGCTGGTGGAGTCTGGGGGAGGCTGAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCATTCTCACTGACTACTGAATGACCT",
                "read frame": 0,
                "strand": True
            },
            "start": 2298801
        },
        "accession": "BN000872",
        "allele": "PG.5.30*01",
        "confirmed": True,
        "end": 2299106,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.5.30",
        "id": "PG.5.30*01-C57BL/6",
        "identified": True,
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "GACCACTCACCATAGACTTTGGGCTCAGCTTTGCTTTCCTTGTTCTTATTTTAAAAGGTAATTCTTAGAAATAAGATCCTGCCTGTTTTGTGTACATGAGAAATAGAAAAATTTGTTTTCTTTCCTCTATTTTGTTTTGTTTTGTTTTGTTAGTGACAGTTTCCAAATCAGCATTCTCTGCTTTGAGGTGCCCAGTGTGAGGTGAAGCTGGTGGAGTCTGGGGGAGGCTGAGTGCAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCATTCTCACTGACTACTGAATGACCT",
            "read frame": 0,
            "strand": True
        },
        "start": 2298802,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "end": 2277088,
            "length": 298,
            "sequence": {
                "nucleotide": "GAGGTGAAGCTGGAGAAGTCTAAGGGGAGGTTTAGTGCAGCCTGGAAGGTCCATGATACTCTACTGTGCAGCCTCTGGATTCACTGTCAGCGACGACTGGTTTGCCTGGGTTTGCCAGGCTCCAAAGAAGGGGCTGGAGATGGGGGGGGGGCATAATTTTTCATGGTGGTGGTAGCCCCTCTTATGCAGACACCTTGAAGAAGTGGGTTGGGCCTAACATATTCAGAATCAATATTTAAAAATTCTAATCCTTGAGTATGTCACTTTTGACCAAGAATATATGAAATATGTTACTGAG",
                "read frame": 0,
                "strand": True
            },
            "start": 2276790
        },
        "accession": "BN000872",
        "allele": "PG.6.32*01",
        "confirmed": True,
        "end": 2277088,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.6.32",
        "id": "PG.6.32*01-C57BL/6",
        "identified": True,
        "length": 297,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "AGGTGAAGCTGGAGAAGTCTAAGGGGAGGTTTAGTGCAGCCTGGAAGGTCCATGATACTCTACTGTGCAGCCTCTGGATTCACTGTCAGCGACGACTGGTTTGCCTGGGTTTGCCAGGCTCCAAAGAAGGGGCTGGAGATGGGGGGGGGGCATAATTTTTCATGGTGGTGGTAGCCCCTCTTATGCAGACACCTTGAAGAAGTGGGTTGGGCCTAACATATTCAGAATCAATATTTAAAAATTCTAATCCTTGAGTATGTCACTTTTGACCAAGAATATATGAAATATGTTACTGAG",
            "read frame": 0,
            "strand": True
        },
        "start": 2276791,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "end": 2102619,
            "length": 305,
            "sequence": {
                "nucleotide": "ATACCACTTGTAGAATCTGAAGGCAGCTTGTTATAGCCCGGAGTGTTCATTAAATTCTCTTCTGAAGCGACTGGATCCACCTTCAGTGATGCACTGGATTCGCCAAGCCCCAGGGAAGGGGCTAGAGTGGGTCAGCAGATATAAAATAATAGTGGAGTGTAAAAACTATGCAGAGTCTGTGAAGGATAGATTTGCCATCTACAGGAACAATTCTAAGAGCTTATTGTAACTACAAATAAACAAATCTAAGAAGTGAGGACACTGCCATGTGTTACTGTGTGAGAGATAAAGTGAGGAAAGTAAAT",
                "read frame": 2,
                "strand": True
            },
            "start": 2102314
        },
        "accession": "BN000872",
        "allele": "PG.8.47*01",
        "confirmed": True,
        "end": 2102619,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.8.47",
        "id": "PG.8.47*01-C57BL/6",
        "identified": True,
        "length": 304,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "TACCACTTGTAGAATCTGAAGGCAGCTTGTTATAGCCCGGAGTGTTCATTAAATTCTCTTCTGAAGCGACTGGATCCACCTTCAGTGATGCACTGGATTCGCCAAGCCCCAGGGAAGGGGCTAGAGTGGGTCAGCAGATATAAAATAATAGTGGAGTGTAAAAACTATGCAGAGTCTGTGAAGGATAGATTTGCCATCTACAGGAACAATTCTAAGAGCTTATTGTAACTACAAATAAACAAATCTAAGAAGTGAGGACACTGCCATGTGTTACTGTGTGAGAGATAAAGTGAGGAAAGTAAAT",
            "read frame": 0,
            "strand": True
        },
        "start": 2102315,
        "strain": "C57BL/6",
        "strand": True
    },
    {
        "BN000872": {
            "end": 2033174,
            "length": 308,
            "sequence": {
                "nucleotide": "GAGATACCACTTGTAGAATCTGAAGGCAGCTTGTTATAGCCTGGAGTGTTCATTAAATTCTCTTCTGAAGCGACTGGATCCACCTTCAGTGATGCACTGGATTCACCAAGCCCCAGGGAAGGGGCTAGAGTGGGTCAGCAGATATAAAATAATAGTGGAGTGTAAAAACTATGCAGAGTCTGTGAAGGATAGATTTGCCATCTACAGGAACAATTCTAAGAGCTTATTGTATCTACAAATAAACAAATCTAAGAAGTGAGGACACTGCCATGTGTTACTGTGTGAGAGATACAGTGAGGAAAGTAAAT",
                "read frame": 0,
                "strand": True
            },
            "start": 2032866
        },
        "accession": "BN000872",
        "allele": "PG.9.52*01",
        "confirmed": True,
        "end": 2033174,
        "family": "unclassified",
        "framed": True,
        "functionality": "P",
        "gene": "PG.9.52",
        "id": "PG.9.52*01-C57BL/6",
        "identified": True,
        "length": 307,
        "organism name": "mus musculus",
        "read frame": 0,
        "region": "VH",
        "sequence": {
            "nucleotide": "AGATACCACTTGTAGAATCTGAAGGCAGCTTGTTATAGCCTGGAGTGTTCATTAAATTCTCTTCTGAAGCGACTGGATCCACCTTCAGTGATGCACTGGATTCACCAAGCCCCAGGGAAGGGGCTAGAGTGGGTCAGCAGATATAAAATAATAGTGGAGTGTAAAAACTATGCAGAGTCTGTGAAGGATAGATTTGCCATCTACAGGAACAATTCTAAGAGCTTATTGTATCTACAAATAAACAAATCTAAGAAGTGAGGACACTGCCATGTGTTACTGTGTGAGAGATACAGTGAGGAAAGTAAAT",
            "read frame": 0,
            "strand": True
        },
        "start": 2032867,
        "strain": "C57BL/6",
        "strand": True
    }
]

notframed = [
    "IGHV5-13-1*01",
    "IGHV(II)-2*01",
    "IGHV(III)-7*01",
    "IGHV5-8-1*01",
    "IGHV6-1*02",
    "IGHV(II)-5*01",
    "IGHV5-11-2*01",
    "IGHV(III)-5*01",
    "IGHV5-12-3*01",
    "IGHV5-8-2*01",
    "IGHV5-7-3*01",
    "IGHV6-2*02",
    "IGHV2-7*02",
    "IGHV5-19*02",
    "IGHV5-7-4*01",
    "IGHV5-21*02",
    "IGHV5-10-1*01",
    "IGHV5-7-2*01",
    "IGHV2-8*02",
    "IGHV5-11-1*01",
    "IGHV5-7-1*01",
    "IGHV5-8-3*01",
    "IGHV6-1-1*01",
    "IGHV(II)-1*01",
    "IGHV(II)-4*01",
    "IGHV15-1*02",
    "IGHV(III)-2*01",
    "IGHV5-18*02",
    "IGHV(III)-3*01",
    "IGHV(III)-4*01",
    "IGHV5-7-5*01",
    "IGHV7-2*02",
    "IGHV5-10-2*01",
    "IGHV(III)-6*01",
    "IGHV(II)-3*01",
    "IGHV5-5-1*01",
    "IGHV5-5*02",
    "IGHV5-13*02",
    "IGHV12-2*02",
    "IGHV2-1*02",
    "IGHV1S57*01",
    "IGHV8S1*01",
    "IGHV8-3*01",
    "IGHV6-2*01",
    "IGHV1-29*01",
    "IGHV1-1*01",
    "IGHV1-17*01",
    "IGHV5-3*01",
    "IGHV1-30*01",
    "IGHV8-15*01",
    "IGHV2-8*01",
    "IGHV1-6*01",
    "IGHV1-3*01",
    "IGHV5-13*01",
    "IGHV6-1*01",
    "IGHV1-62*01",
    "IGHV1-45*01",
    "IGHV1-44*01",
    "IGHV10-4*01",
    "IGHV5-7*01",
    "IGHV5-10*01",
    "IGHV(III)-8*01",
    "IGHV5-19*01",
    "IGHV1-68*01",
    "IGHV1-38*01",
    "IGHV3-2*01",
    "IGHV2-1*01",
    "IGHV8-16*01",
    "IGHV10-2*01",
    "IGHV1-2*01",
    "IGHV(III)-1*01",
    "IGHV(III)-10*01",
    "IGHV5-5*01",
    "IGHV1-33*01",
    "IGHV5-18*01",
    "IGHV1-40*01",
    "IGHV(III)-9*01",
    "IGHV8-1*01",
    "IGHV1-57*01",
    "IGHV5-8*01",
    "IGHV1-65*01",
    "IGHV5-11*01",
    "IGHV1-73*01",
    "IGHV1-10*01",
    "IGHV1-41*01",
    "IGHV15-1*01",
    "IGHV1-86*01",
    "IGHV(III)-11*01"
]

oops = {
    "accession": "BN000872",
    "allele": "IGHV6-2*01",
    "confirmed": True,
    "end": 1981725,
    "family": "J606",
    "framed": False,
    "functionality": "P",
    "gene": "IGHV6-2",
    "id": "IGHV6-2*01-C57BL/6",
    "identified": True,
    "imgt allele name": "IGHV6-2*01",
    "imgt family name": "IGHV6",
    "imgt gene name": "IGHV6-2",
    "length": 295,
    "organism name": "mus musculus",
    "read frame": 0,
    "region": "VH",
    "sequence": {
        "nucleotide": "GAGTTGCAGCTTATGGAGTCTTGGGGAAGCTTGTTAAAGCTCCAGGGTTCTGTAAGACATTCTTGTGCAGCCGCTGGATTCACTTTCAGAGAATACTATATCAGATGGGTCCTGTAGATTTTAGAAAAAGATCTTGAGTCCTTGATTTAATTACAAACACAGCTGGGGGGTGATTACAGAGTATGCTTCCTATGTGAGAAGGGCACTTCCCATTTCAAGAGATCATACAAAAAATAAAATTTCATGCCTCCAAATAAACAGATTGGAGACTTCTTTTTTCCCTGAGAAGAAATGG",
        "read frame": 0,
        "strand": True
    },
    "start": 1981430,
    "strain": "C57BL/6",
    "strand": True,
    "verified": True
}

for gene in good:
    index['good'][gene['id']] = gene

for gene in bad:
    index['bad'][gene['id']] = gene


for k,v in index['good'].items():
    b = index['bad'][k]
    g = v

    g['length'] = g['end'] - g['start']
    if 'framed' in b: g['framed'] = b['framed']
    if not g['framed']:
        if 'read frame' in g:
            del g['read frame']
        g['sequence']['read frame'] = 0

    if 'BN000872' in b:
        g['BN000872'] = b['BN000872']
        g['BN000872']['length'] = g['BN000872']['end'] - g['BN000872']['start']
        if not g['framed']:
            g['BN000872']['read frame'] = 0
        else:
            gap = g['BN000872']['start'] - g['start']
            g['BN000872']['sequence']['read frame'] = (3 - (gap - g['read frame'])%3)%3
    if 'strand' not in g:
        print(k)


# print(to_json(good))
