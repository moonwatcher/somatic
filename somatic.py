#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import logging
import re
import json
import io

from subprocess import Popen, PIPE
from io import StringIO, BytesIO
from pymongo import MongoClient, DESCENDING, ASCENDING
from argparse import ArgumentParser
from bson.objectid import ObjectId
from datetime import timedelta, datetime

def to_json(node):
    def handler(o):
        result = None
        if isinstance(o, datetime):
            result = o.isoformat()
        if isinstance(o, ObjectId):
            result = str(o)
        if isinstance(o, set):
            result = list(o)
        if isinstance(o, Sequence):
            result = o.node
        return result

    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4, default=handler)

log_levels = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

expression = {
    'complement': { 
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'N': 'N' 
    }, 
    'nucleic to amino': {
        'GCT':'A',
        'GCC':'A',
        'GCA':'A',
        'GCG':'A',
        'TTA':'L',
        'TTG':'L',
        'CTT':'L',
        'CTC':'L',
        'CTA':'L',
        'CTG':'L',
        'CGT':'R',
        'CGC':'R',
        'CGA':'R',
        'CGG':'R',
        'AGA':'R',
        'AGG':'R',
        'AAA':'K',
        'AAG':'K',
        'AAT':'N',
        'AAC':'N',
        'ATG':'M',
        'GAT':'D',
        'GAC':'D',
        'TTT':'F',
        'TTC':'F',
        'TGT':'C',
        'TGC':'C',
        'CCT':'P',
        'CCC':'P',
        'CCA':'P',
        'CCG':'P',
        'CAA':'Q',
        'CAG':'Q',
        'TCT':'S',
        'TCC':'S',
        'TCA':'S',
        'TCG':'S',
        'AGT':'S',
        'AGC':'S',
        'GAA':'E',
        'GAG':'E',
        'ACT':'T',
        'ACC':'T',
        'ACA':'T',
        'ACG':'T',
        'GGT':'G',
        'GGC':'G',
        'GGA':'G',
        'GGG':'G',
        'TGG':'W',
        'CAT':'H',
        'CAC':'H',
        'TAT':'Y',
        'TAC':'Y',
        'ATT':'I',
        'ATC':'I',
        'ATA':'I',
        'GTT':'V',
        'GTC':'V',
        'GTA':'V',
        'GTG':'V',
        'TAA':'*',
        'TGA':'*',
        'TAG':'*',
    },
    'buffer size': 1,
    'regions': ['VH', 'DH', 'JH'],
    'fasta line length': 80,
    'igblast result start': re.compile('^# IGBLASTN '),
    'igblast query': re.compile('^# Query: (P:<query_id>.*)$'),
    'igblast hit table': re.compile('^# Hit table'),
    'igblast result end': re.compile('^# BLAST processed '),
    'igblast reversed query': '# Note that your query represents the minus strand of a V gene and has been converted to the plus strand.',
    'nucleotide sequence': re.compile('^[atcgnATCGN]+$'),
    'minimum v identity': 0.7, 
    'minimum d identity': 0.7, 
    'minimum j identity': 0.7, 
    'minimum v alignment': 40, 
    'minimum d alignment': 4, 
    'minimum j alignment': 30, 
    'd overlap penalty factor': 0.5, 
    'igblast hit': re.compile(
        r"""
        (?P<region>[VDJ])\t
        (?:(?P<query_polarity>reversed)\|)?(?:[^,]+)\t
        (?P<subject_id>[^,]+)\t
        (?P<query_start>[^,]+)\t
        (?P<query_end>[^,]+)\t
        (?P<subject_start>[^,]+)\t
        (?P<subject_end>[^,]+)\t
        (?P<gap_openings>[^,]+)\t
        (?P<gaps>[^,]+)\t
        (?P<mismatch>[^,]+)\t
        (?P<identical>[^,]+)\t
        (?P<bit_score>[^,]+)\t
        (?P<evalue>[^,]+)\t
        (?P<alignment_length>[^,]+)\t
        (?P<subject_strand>[^,]+)
        """,
        re.VERBOSE
    ),
    'imgt fasta header': re.compile(
        r"""
        ^>
        (?P<imgt_accession>[^|]*)\|
        (?P<allele_name>(?P<gene_name>(?P<subgroup>[^-]+)[^\*]+)\*[0-9]+)\|
        (?P<organism_name>[^|_]*)(?:_(?P<strain>[^|]+))?\|
        (?P<functionality>[\[\(]?(?:F|ORF|P)[\]\)]?)\|
        (?P<region>V|D|J)-REGION\|
        (?:(?P<start>[0-9]+)\.\.(?P<end>[0-9]+))?\s*\|
        (?P<length>[0-9]+)\snt\|
        (?P<read_frame>[123]|NR)\|
        [^|]*\|
        [^|]*\|
        [^|]*\|
        [^|]*\|
        [^|]*\|
        [^|]*\|
        (?P<polarity>rev-compl)?\s*\|$
        """,
        re.VERBOSE
    ),
    'expand hit': re.compile(
        r"""
        (?P<region>VH|DH|JH),
        (?P<subject_id>[^,]+),
        (?P<query_start>[^,]+),
        (?P<query_end>[^,]+),
        (?P<subject_start>[^,]+),
        (?P<subject_end>[^,]+),
        (?P<gap_openings>[^,]+),
        (?P<gaps>[^,]+),
        (?P<mismatch>[^,]+),
        (?P<identical>[^,]+),
        (?P<bit_score>[^,]+),
        (?P<evalue>[^,]+),
        (?P<alignment_length>[^,]+),
        (?P<subject_strand>[^,]),
        (?P<query_strand>[^,])
        """,
        re.VERBOSE
    ),
    'igblast command': [
        'igblastn',
        '-germline_db_V', 'database/mouse_imgt_vh',
        '-germline_db_J', 'database/mouse_imgt_jh',
        '-germline_db_D', 'database/mouse_imgt_dh',
        '-organism', 'mouse',
        '-domain_system', 'imgt',
        '-query', '-',
        '-auxiliary_data', 'optional_file/mouse_gl.aux',
        '-show_translation',
        '-outfmt',
        '7 qseqid sseqid qstart qend sstart send gapopen gaps mismatch pident bitscore evalue length sstrand',
    ],
    'igblast compressed hit': '{region},{subject id},{query start},{query end},{subject start},{subject end},{gap openings},{gaps},{mismatch},{identical},{bit score},{evalue},{alignment length},{subject strand},{query strand}',
}


class Sequence(object):
    def __init__(self, node=None):
        self.log = logging.getLogger('Sequence')
        self._reversed = None
        if node is None:
            self.node = {
                'read frame': 0,
                'strand': True
            }
        else:
            self.node = node

    def __str__(self):
        buffer = StringIO()
        if self.read_frame > 0:
            buffer.write(self.nucleotide[0:self.read_frame])
            buffer.write(' ')
        for index in range(self.read_frame, self.length, 3):
            buffer.write(self.nucleotide[index:index + 3])
            buffer.write(' ')
        buffer.write('\n')
        buffer.seek(0)
        return buffer.read()

    def reset(self):
        for k in ['codon']:
            if k in self.node: del self.node[k]
        self._reversed = None

    def clone(self):
        clone = {
            'nucleotide': self.nucleotide, 
            'read frame': self.read_frame,
            'strand': self.strand
        }
        if 'quality' in self.node:
            clone['quality'] = self.quality
        return Sequence(clone)

    def crop(self, start, end):
        cropped = {
            'nucleotide': self.nucleotide[start:end], 
            'read frame': (3 - (start - self.read_frame)%3)%3,
            'strand': self.strand
        }
        if 'quality' in self.node:
            cropped['quality'] = self.quality[start:end]
        return Sequence(cropped)

    @property
    def read_frame(self):
        return self.node['read frame']

    @read_frame.setter
    def read_frame(self, value):
        self.node['read frame'] = value
        self.reset()

    @property
    def strand(self):
        return self.node['strand']

    @strand.setter
    def strand(self, value):
        self.node['strand'] = value
        self.reset()

    @property
    def nucleotide(self):
        if 'nucleotide' in self.node:
            return self.node['nucleotide']
        else:
            return None

    @nucleotide.setter
    def nucleotide(self, value):
        self.node['nucleotide'] = value
        if not self.node['nucleotide']:
            del self.node['nucleotide']
        self.reset()

    @property
    def quality(self):
        if 'quality' in self.node:
            return self.node['quality']
        else:
            return None

    @quality.setter
    def quality(self, value):
        self.node['quality'] = value
        if not self.node['quality']:
            del self.node['quality']
        self.reset()

    @property
    def codon(self):
        if 'codon' not in self.node and self.nucleotide:
            start = self.read_frame
            end = start + 3
            codon = []
            while(not end > self.length):
                acid = self.nucleotide[start:end]
                if 'N' in acid:
                    codon.append('X')
                else:
                    codon.append(expression['nucleic to amino'][acid])
                start = end
                end = start + 3
            if codon:
                self.node['codon'] = ''.join(codon)
        if 'codon' in self.node:
            return self.node['codon']
        else:
            return None

    @property
    def length(self):
        if self.nucleotide:
            return len(self.nucleotide)
        else:
            return None

    @property
    def reversed(self):
        if self._reversed is None and self.nucleotide:
            reversed = {
                'nucleotide': ''.join([ expression['complement'][b] for b in list(self.nucleotide) ][::-1]), 
                'read frame': (self.length - self.read_frame) % 3,
                'strand': not self.strand
            }
            if 'quality' in self.node:
                reversed['quality'] = self.quality[::-1]
                
            self._reversed = Sequence(reversed)
        return self._reversed


class Reference(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Reference')
        self.pipeline = pipeline
        self.node = node
        
        if self.node is None: self.node = {}
            
        if 'sequence' in self.node:
            self.sequence = Sequence(self.node['sequence'])
        else:
            self.sequence = Sequence()
            self.node['sequence'] = self.sequence.node

    @property
    def id(self):
        return self.node['allele name']

    @property
    def region(self):
        return self.node['region']

    @property
    def strain(self):
        if 'strain' in self.node:
            return self.node['strain']
        else:
            return None

    @property
    def framed(self):
        return self.node['framed']

    @property
    def functionality(self):
        result = None
        if 'functionality' in self.node:
            if 'F' in self.node['functionality']:
                result = 'F'
            if 'P' in self.node['functionality']:
                 result = 'P'
            if 'ORF' in self.node['functionality']:
                result = 'O'
        return result


class Sample(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Sample')
        self.pipeline = pipeline
        self.node = node
        
        if self.node is None: self.node = {}
            
        if 'sequence' in self.node:
            self.sequence = Sequence(self.node['sequence'])
        else:
            self.sequence = Sequence()
            self.node['sequence'] = self.sequence.node

    @property
    def id(self):
        return self.node['id']

    @id.setter
    def id(self, value):
        self.node['id'] = value

    @property
    def library(self):
        return self.node['library']

    @library.setter
    def library(self, value):
        self.node['library'] = value

    @property
    def valid(self):
        return self.node['valid']

    @valid.setter
    def valid(self, value):
        self.node['valid'] = value

    def reverse(self):
        if self.sequence:
            self.sequence = self.sequence.reversed
            self.node['sequence'] = self.sequence.node

    @property
    def hit(self):
        if 'hit' not in self.node:
            self.node['hit'] = []
        return self.node['hit']

    @property
    def region(self):
        if 'region' not in self.node:
            self.node['region'] = {}
        return self.node['region']

    def expand(self):
        if 'hit' in self.node:
            for hit in self.node['hit']:
                if 'compressed hit' in hit:
                    match = expression['expand hit'].search(hit['compressed hit'])
                    if match:
                        parsed = dict(((k.replace('_', ' '),v) for k,v in match.groupdict().items()))
                        for k,v in parsed.items():
                            if k in [
                                'query start',
                                'query end',
                                'subject start',
                                'subject end',
                                'gap openings',
                                'gaps',
                                'mismatch',
                                'alignment length'
                            ]:
                                try:
                                    hit[k] = int(parsed[k])
                                except ValueError as e:
                                    self.log.warning('could not parse value %s for %s as int', parsed[k], k)
                            elif k in [
                                'identical',
                                'bit score',
                                'evalue',
                            ]:
                                try:
                                    hit[k] = float(parsed[k])
                                except ValueError:
                                    self.log.warning('could not parse value %s for %s as float', parsed[k], k)
                            else:
                                hit[k] = v

    def compress(self):
        if 'hit' in self.node:
            for hit in self.node['hit']:
                try:
                    hit['compressed hit'] = expression['igblast compressed hit'].format(**hit)
                except KeyError as e:
                    self.log.error('could not compress hit because %s was missing', e)

    @property
    def gapped(self):
        if 'gapped' not in self.node:
            self.node['gapped'] = False
            for hit in self.hit:
                if hit['gap openings'] > 0:
                    self.node['gapped'] = True
                    break
        return self.node['gapped']

    def _complement_from_reference(self):
        for hit in self.hit:
            hit['picked'] = False
            reference = self.pipeline.reference_for(hit['subject id'])
            hit['framed'] = reference.framed
            hit['functionality'] = reference.functionality
            if reference.strain:
                hit['strain'] = reference.strain
            hit['subject'] = reference.sequence.crop(hit['subject start'], hit['subject end'])
            hit['query'] = self.sequence.crop(hit['query start'], hit['query end'])
        return True

    def _pick_jh_region(self):
        picked = False
        for hit in self.hit:
            if (hit['framed'] and 
                hit['region'] == 'JH' and 
                hit['identical'] >=  expression['minimum j identity'] and
                hit['alignment length'] >= expression['minimum j alignment']):
                
                self.region['JH'] = hit
                hit['picked'] = True
                self.sequence.read_frame = (hit['query start'] + hit['subject'].read_frame) % 3
                picked = True
                break
        return picked

    def _assign_frame(self):
        for hit in self.hit:
            hit['query'].read_frame = 2 - (hit['query start'] - self.sequence.read_frame - 1) % 3
            if hit['framed']:
                if hit['query'].read_frame == hit['subject'].read_frame:
                    hit['in frame'] = 'Y'
                else:
                    hit['in frame'] = 'N'
            else:
                hit['in frame'] = '*'
        return True

    def _pick_vh_region(self):
        picked = False
        for hit in self.hit:
            if (hit['framed'] and 
                hit['region'] == 'VH' and 
                hit['identical'] >=  expression['minimum v identity'] and
                hit['alignment length'] >= expression['minimum v alignment']):
                
                self.region['VH'] = hit
                hit['picked'] = True
                picked = True
                break
        return picked

    def _pick_dh_region(self):
        picked = False
        possible = []
        for hit in self.hit:
            if hit['region'] == 'DH':
                h = {}
                for k,v in hit.items():
                    if k not in ('query', 'subject'): h[k] = v
                reference = self.pipeline.reference_for(hit['subject id'])
                    
                h['overlap'] = 0
                if self.region['VH']['query end'] > hit['query start']:
                    h['query start'] = self.region['VH']['query end']
                    overlap = h['query start'] - hit['query start']
                    h['subject start'] += overlap
                    h['overlap'] += overlap
                    
                if hit['query end'] > self.region['JH']['query start']:
                    h['query end'] = self.region['JH']['query start']
                    overlap = hit['query end'] - h['query end']
                    h['subject end'] -= overlap
                    h['overlap'] += overlap
                    
                h['score'] = hit['bit score'] - h['overlap'] * expression['d overlap penalty factor']
                h['alignment length'] = h['query end'] - h['query start']
                h['query'] = self.sequence.crop(h['query start'], h['query end'])
                h['subject'] = reference.sequence.crop(h['subject start'], h['subject end'])
                
                if h['query end'] > h['query start']:
                    similar = 0
                    different = 0
                    longest = 0
                    stretch = 0
                    for index in range(h['query'].length):
                        if h['subject'].nucleotide[index] == h['query'].nucleotide[index]:
                            similar += 1
                            stretch += 1
                            longest = max(longest, stretch)
                        else:
                            different += 1
                            stretch = 0
                        
                    h['identical'] = float(similar) / float(h['alignment length'])
                    if (longest >= expression['minimum d alignment'] and
                        h['identical'] >= expression['minimum d identity']):
                        possible.append(h)

        if possible:
            possible.sort(reverse=True, key=lambda x: x['score'])
            self.region['DH'] = possible[0]
            for h in self.hit:
                if h['region'] == 'DH' and h['subject id'] == self.region['DH']['subject id']:
                    h['picked'] = True
            picked = True
        return picked

    def _identify_cdr3(self):
        start = None
        end = None
        
        offset = None
        # look for the most upstream cycteine on the VH region
        for index,codon in enumerate(reversed(self.region['VH']['query'].codon)):
            if codon == 'C':
                offset = self.region['VH']['query'].read_frame + (len(self.region['VH']['query'].codon) - index - 1) * 3
                break
                
        if offset is not None:
            start = self.region['VH']['query start'] + offset
            
        offset = None
        # look for the most downstream tryptophan JH region
        for index,codon in enumerate(self.region['JH']['query'].codon):
            if codon == 'W':
                offset = self.region['JH']['query'].read_frame + (index + 1) * 3
                break
                
        if offset is not None:
            end = self.region['JH']['query start'] + offset
        
        if start and end:
            sequence = Sequence({
                'nucleotide': self.sequence.nucleotide[start:end],
                'strand': self.sequence.strand,
                'read frame': 0
            })
            sequence.codon
            self.region['CDR3'] = {
                'query start': start,
                'query end': end,
                'region': 'CDR3',
                'query': sequence
            }
            

    def _check_for_stop_codon(self):
        if '*' in self.sequence.codon:
            self.node['premature termination'] = True
        else:
            self.node['premature termination'] = False
        return not self.node['premature termination']

    def _check_for_aligned_frames(self):
        if self.region['JH'] and self.region['VH']:
            if self.region['VH']['in frame'] == 'Y' and self.region['JH']['in frame'] == 'Y':
                self.node['in frame'] = True
            else:
                self.node['in frame'] = False
        else:
            self.node['in frame'] = False
        return True

    def analyze(self):
        valid = not self.gapped
        valid = valid and self._complement_from_reference()
        valid = valid and self._pick_jh_region()
        valid = valid and self._assign_frame()
        valid = valid and self._pick_vh_region()
        valid = valid and self._check_for_aligned_frames()
        valid and self._pick_dh_region()
        valid and self._check_for_stop_codon()
        valid and self._identify_cdr3()
        self.valid = valid

    def to_json(self):
        def handler(o):
            result = None
            if isinstance(o, datetime):
                result = o.isoformat()
            if isinstance(o, ObjectId):
                result = str(o)
            if isinstance(o, set):
                result = list(o)
            if isinstance(o, Sequence):
                result = o.node
            return result
            
        return json.dumps(self.node, sort_keys=False, ensure_ascii=False, indent=4, default=handler)

    def to_summary(self):
        if self.hit:
            buffer = StringIO()
            buffer.write(self.id)
            buffer.write(' | ')
            for r in ('VH', 'DH', 'JH'):
                if r in self.region:
                    buffer.write(r)
                    buffer.write(' : ')
                    buffer.write(self.region[r]['subject id'])
                    buffer.write(' | ')
            
            if self.node['in frame']:
                buffer.write('in frame')
            else:
                buffer.write('out of frame')
                
            if self.node['premature termination']:
                buffer.write(' | prematurely terminated')
                
            buffer.seek(0)
            return buffer.read()
        else:
            return ''

    def to_alignment_diagram(self):
        if self.hit:
            buffer = StringIO()
            subject_length = max(max([ len(hit['subject id']) for hit in self.hit ]), len('subject'))
            alignment_offset = subject_length + 12
            
            # print a title
            buffer.write('{: <{}}{}'.format('Summary', alignment_offset, self.to_summary()))
            buffer.write('\n')
            
            # print the coordinate system
            buffer.write('{: <{}} '.format('Subject', subject_length))
            buffer.write('{: <{}} '.format('R', 4))
            buffer.write('{: <{}} '.format('F', 1))
            buffer.write('{: <{}} '.format('I', 1))
            buffer.write('{: <{}} '.format('P', 1))
            if self.sequence.read_frame > 0:
                buffer.write('{: <{}} '.format(0, self.sequence.read_frame))
            for index in range(self.sequence.read_frame, self.sequence.length, 3):
                buffer.write('{: <3} '.format(index))
            buffer.write('\n')
            
            # print the sample nucleotide sequence
            buffer.write(' ' * alignment_offset)
            if self.sequence.read_frame > 0:
                buffer.write(self.sequence.nucleotide[0:self.sequence.read_frame])
                buffer.write(' ')
            for index in range(self.sequence.read_frame, self.sequence.length, 3):
                buffer.write(self.sequence.nucleotide[index:index + 3])
                buffer.write(' ')
            buffer.write('\n')
            
            # print the sample quality sequence
            buffer.write(' ' * alignment_offset)
            if self.sequence.read_frame > 0:
                buffer.write(self.sequence.quality[0:self.sequence.read_frame])
                buffer.write(' ')
            for index in range(self.sequence.read_frame, self.sequence.length, 3):
                buffer.write(self.sequence.quality[index:index + 3])
                buffer.write(' ')
            buffer.write('\n')
            
            # print the sample codon sequence
            buffer.write(' ' * alignment_offset)
            if self.sequence.read_frame > 0:
                buffer.write('{: <{}} '.format(' ', self.sequence.read_frame))
            for codon in self.sequence.codon:
                buffer.write(codon)
                buffer.write('   ')
            buffer.write('\n')
            
            for hit in self.hit:
                buffer.write('{: <{}} '.format(hit['subject id'], subject_length))
                buffer.write('{: <{}} '.format(hit['region'], 4))
                buffer.write('{: <{}} '.format(hit['functionality'], 1))
                buffer.write('{: <{}} '.format(hit['in frame'] if 'in frame' in hit else ' ', 1))
                buffer.write('{: <{}} '.format('Y' if hit['picked'] else 'N', 1))
                
                subject = hit['subject'].clone()
                subject.read_frame = hit['query'].read_frame
                
                hit_offset = 0
                if hit['query start'] > 0:
                    if self.sequence.read_frame > 0:
                        hit_offset = (int((hit['query start'] - self.sequence.read_frame) / 3)) + hit['query start'] + 1
                    else:
                        hit_offset = (int(hit['query start'] / 3)) + hit['query start']
                    buffer.write(' ' * hit_offset)
                    
                # nucleotide mask
                try:
                    mask = []
                    for index,n in enumerate(subject.nucleotide):
                        if n == hit['query'].nucleotide[index]: mask.append('-')
                        else: mask.append(n)
                    nucleotide_mask = ''.join(mask)
                except IndexError as e:
                    print(e)
                    print(to_json(hit))
                    print(to_json(self.node))
                    
                # print nucleotides
                if subject.read_frame > 0:
                    buffer.write(nucleotide_mask[0:subject.read_frame])
                    buffer.write(' ')
                    
                for index in range(subject.read_frame, subject.length, 3):
                    buffer.write(nucleotide_mask[index: index + 3])
                    buffer.write(' ')
                buffer.write('\n')
                
                # print codons sequence if anything is different
                display = False
                offset = int((hit['query start'] - self.sequence.read_frame + subject.read_frame) / 3)
                mask = []
                for index,codon in enumerate(subject.codon):
                    if codon == hit['query'].codon[index]:
                        mask.append('-')
                    else:
                        mask.append(codon)
                        display = True
                codon_mask = ''.join(mask)
                
                if display:
                    buffer.write(' ' * alignment_offset)
                    if hit['query start'] > 0:
                        buffer.write(' ' * hit_offset)
                        
                    if subject.read_frame > 0:
                        buffer.write('{: <{}} '.format(' ', subject.read_frame))
                        
                    for codon in mask:
                        buffer.write(codon)
                        buffer.write('   ')
                    buffer.write('\n')
            buffer.seek(0)
            return buffer.read()
        else:
            return ''

    def view(self):
        if self.valid:
            print(self.to_alignment_diagram())

    def info(self):
        print(self.to_json())

    def __str__(self):
        return self.to_alignment_diagram()


class Block(object):
    def __init__(self, pipeline):
        self.log = logging.getLogger('Block')
        self.pipeline = pipeline
        self.reset()

    def __str__(self):
        buffer = StringIO()
        for sample in self.buffer:
            if sample.valid:
                buffer.write(str(sample))
                buffer.write('\n')
        buffer.seek(0)
        return buffer.read()

    @property
    def buffer(self):
        return self.node['buffer']

    @property
    def document(self):
        return self.node['document']

    @property
    def lookup(self):
        return self.node['lookup']

    @property
    def size(self):
        return len(self.buffer)

    @property
    def empty(self):
        return self.size == 0

    @property
    def full(self):
        return not self.size < expression['buffer size']

    def reset(self):
        self.node = {
            'buffer': [],
            'document': [],
            'lookup': {}
        }

    def add(self, sample):
        self.buffer.append(sample)
        self.lookup[sample.id] = sample
        self.document.append(sample.node)

    def fill(self, library):
        if self.read():
            buffer = self.search()
            if buffer:
                self.parse_igblast(buffer)
                for sample in self.buffer:
                    sample.library = library
                    sample.analyze()
        return not self.empty

    def read(self):
        self.reset()
        state = None
        sample = Sample(self.pipeline)
        for line in sys.stdin:
            if line:
                line = line.strip()
                if state == None:
                    if line[0] == '@':
                        sample.id = line
                        state = 'id'
                        
                elif state == 'id':
                    sample.sequence.nucleotide = line
                    state = 'sequence'
                    
                elif state == 'sequence':
                    if line[0] == '+':
                        state = 'redundant'
                        
                elif state == 'redundant':
                    sample.sequence.quality = line
                    self.add(sample)
                    state = None
                    if self.full:
                        break
                    else:
                        sample = Sample(self.pipeline)
            else:
                break
                
        return not self.empty

    def search(self):
        result = None
        process = Popen(
            args=expression['igblast command'],
            cwd='/Users/lg/code/somatic/igblast',
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate(input=self.to_fasta().read().encode('utf8'))
        if output:
            result = StringIO(output.decode('utf8'))
        return result

    def to_fasta(self):
        buffer = StringIO()
        for sample in self.buffer:
            buffer.write('>{}\n'.format(sample.id))
            begin = 0
            end = 0
            while end < sample.sequence.length:
                end = min(begin + expression['fasta line length'], sample.sequence.length)
                buffer.write(sample.sequence.nucleotide[begin:end])
                buffer.write('\n')
                begin = end
        buffer.seek(0)
        return buffer

    def parse_igblast(self, buffer):
        state = 0
        sample = None
        for line in buffer:
            line = line.strip()
            if state == 0 and line.startswith('# IGBLASTN '):
                state = 1
            elif state == 1 and line.startswith('# Query: '):
                id = line[9:]
                if id in self.lookup:
                    sample = self.lookup[id]
                    state = 2
                else:
                    self.log.error('could not locate %s in buffer', id)
                    state = 0
            elif state == 2 and line.startswith('# Hit table '):
                state = 3
            elif state == 2 and line.startswith(expression['igblast reversed query']):
                sample.sequence.strand = False
                sample.reverse()
            elif (state == 3 or state == 4) and not line.startswith('#'):
                state = 4
                match = expression['igblast hit'].search(line)
                if match:
                    hit = dict(((k.replace('_', ' '),v) for k,v in match.groupdict().items()))
                    if 'query polarity' in hit and hit['query polarity'] == 'reversed':
                        if hit['subject strand'] == 'plus':
                            hit['query strand'] = 'minus'
                        else:
                            hit['query strand'] = 'plus'
                        del hit['query polarity']
                    else:
                        hit['query strand'] = hit['subject strand']
                        
                    for k in [
                        'query strand',
                        'subject strand'  
                    ]:
                        if hit[k] == 'plus': hit[k] = True
                        elif hit[k] == 'minus': hit[k] = False
                        else:
                            self.log.warning('could not parse value %s for %s as a strand', hit[k], k)
                            hit[k] = None
                            
                    for k in [
                        'identical',
                        'bit score',
                        'evalue',
                    ]: 
                        try: hit[k] = float(hit[k])
                        except ValueError:
                            self.log.warning('could not parse value %s for %s as int', hit[k], k)
                            hit[k] = None
                        
                    for k in [
                        'subject end',
                        'gap openings',
                        'subject start',
                        'mismatch',
                        'query start',
                        'alignment length',
                        'gaps',
                        'query end',
                    ]:
                        try: hit[k] = int(hit[k])
                        except ValueError as e:
                            self.log.warning('could not parse value %s for %s as int', hit[k], k)
                            hit[k] = None

                    # fix the region for heavy chain
                    if 'region' in hit: hit['region'] += 'H'
                    
                    if 'identical' in hit: hit['identical'] /= 100.0
                    if 'subject start' in hit: hit['subject start'] -= 1
                    if 'query start' in hit: hit['query start'] -= 1
                    if 'read frame' in hit: hit['read frame'] -= 1
                    sample.hit.append(hit)
                    
            elif state == 4 and line.startswith('#'):
                state = 0
                sample = None
                if line.startswith('# IGBLASTN '):
                    state = 1

    def view(self):
        for sample in self.buffer:
            sample.view()

    def info(self):
        for sample in self.buffer:
            sample.info()

    def simulate(self, json, alignment):
        for sample in self.buffer:
            if json: sample.info()
            if alignment: sample.view()


class Pipeline(object):
    def __init__(self):
        self.log = logging.getLogger('Pipeline')
        self.count = 0
        self._connection = None
        self._reference_sequence = {}

    @property
    def connection(self):
        if self._connection is None:
            try:
                self._connection = MongoClient('mongodb://localhost')
            except pymongo.errors.ConnectionFailure as e:
                self.log.error('failed to establish connection %s', e)
            else:
                self.log.debug('connection established')
        return self._connection

    @property
    def database(self):
        if self.connection is not None:
            return self.connection['somatic']
        else:
            return None

    def close(self):
        if self._connection is not None:
            self._connection.close()

    def reference_for(self, id):
        if id not in self._reference_sequence:
            node = self.database['reference'].find_one({'allele name': id})
            if node:
                reference = Reference(self, node)
                self._reference_sequence[id] = reference
            else:
                self._reference_sequence[id] = None
        return self._reference_sequence[id]

    def rebuild(self):
        index = [
            {
                'collection': 'reference',
                'index': [
                    { 'key': [( 'allele name', ASCENDING )], 'unique': True, 'name': 'reference allele name' },
                    { 'key': [( 'region', ASCENDING )], 'unique': False, 'name': 'reference region' },
                    { 'key': [( 'subgroup', ASCENDING )], 'unique': False, 'name': 'reference subgroup' },
                    { 'key': [( 'gene name', ASCENDING )], 'unique': False, 'name': 'reference gene name' }
                ]
            }
        ]
        for table in index:
            collection = self.database[table['collection']]
            existing = collection.index_information()
            for definition in table['index']:
                print(definition)
                if definition['name'] in existing:
                    self.log.info('dropping index %s on collection %s', definition['name'], table['collection'])
                    collection.drop_index(definition['name'])
                    
                self.log.info('rebuilding index %s on collection %s', definition['name'], table['collection'])
                collection.create_index(definition['key'], name=definition['name'], unique=definition['unique'])

    def populate_reference(self, path):
        def save_reference_sample(collection, sample):
            if sample is not None and sample['sequence']:
                sample['sequence'] = {
                    'nucleotide': sample['sequence'],
                    'strand': sample['strand'],
                    'read frame': 0
                }
                del sample['strand']
                if 'read frame' in sample:
                    sample['sequence']['read frame'] = sample['read frame']
                    del sample['read frame']
                collection.save(sample)
                
        collection = self.database['reference']
        sample = None
        with io.open(path, 'rb') as fasta:
            for line in fasta:
                if line:
                    line = line.strip().decode('utf8')
                    if line[0] == '>':
                        # at the start of a new record save the completed one
                        save_reference_sample(collection, sample)
                            
                        # initialize a new sample
                        sample = { 'id': line, 'sequence': '' }
                        match = expression['imgt fasta header'].search(line)
                        if match:
                            parsed = dict(((k.replace('_', ' '),v) for k,v in match.groupdict().items()))
                            for k,v in parsed.items():
                                if v:
                                    try:
                                        if k in ['start', 'end', 'length']:
                                            v = int(v)
                                        elif k == 'read frame':
                                            if v == 'NR':
                                                v = None
                                            else:
                                                v = int(v)
                                    except ValueError as e:
                                        self.log.error('unable to decode %s as int for %s', v, k)
                                    else:
                                        if v is not None:
                                            sample[k] = v
                                            
                            # fix the region for heavy chain
                            sample['region'] += 'H'
                            
                            # fix the strand
                            sample['strand'] = True
                            if 'polarity' in sample:
                                if sample['polarity'] == 'rev-compl':
                                    sample['strand'] = False
                                    del sample['polarity']
                                    
                            # fix the start offset
                            if 'start' in sample: sample['start'] -= 1
                            
                            # fix the read frame
                            if 'read frame' in sample:
                                sample['read frame'] -= 1 
                                sample['framed'] = True
                            else:
                                sample['framed'] = False
                    else:
                        if expression['nucleotide sequence'].search(line):
                            sample['sequence'] += line.upper()
                else:
                    break
                    
            # save the last one if it has a sequence
            save_reference_sample(collection, sample)

    def reference_to_blast_fasta(self, region):
        buffer = StringIO()
        collection = self.database['reference']
        query = { 'region': region, 'sequence.strand': True }
        cursor = collection.find(query)
        for sample in cursor:
            buffer.write('>')
            buffer.write(sample['allele name'])
            buffer.write('\n')
            begin = 0
            end = 0
            while end < sample['length']:
                end = min(begin + expression['fasta line length'], sample['length'])
                buffer.write(sample['sequence']['nucleotide'][begin:end])
                buffer.write('\n')
                begin = end
        buffer.seek(0)
        print(buffer.read())

    def reference_to_auxiliary(self, region):
        buffer = StringIO()
        collection = self.database['reference']
        query = { 'region': region, 'sequence.strand': True }
        cursor = collection.find(query)
        for sample in cursor:
            if sample['framed']:
                buffer.write(sample['allele name'])
                buffer.write('\t')
                buffer.write(str(sample['sequence']['read frame']))
                buffer.write('\t')
                buffer.write(sample['region'])
                buffer.write('\n')
        buffer.seek(0)
        print(buffer.read())

    def populate(self, library):
        block = Block(self)
        collection = self.database[library]
        while block.fill(library):
            result = collection.insert_many(block.document)
            self.log.debug('%s so far', self.count)
            self.count += block.size

    def view(self, query):
        collection = self.database['sample']
        cursor = collection.find(query)
        for sample in cursor:
            sample.view()

    def info(self, query):
        collection = self.database['sample']
        cursor = collection.find(query)
        for sample in cursor:
            sample.info()

    def simulate(self, library, json, alignment):
        block = Block(self)
        while block.fill(library):
            block.simulate(json, alignment)


def decode_cli():
    env = {}
    
    # -- global arguments for all actions --
    p = ArgumentParser()
    p.add_argument(
        '-v', '--verbosity',
        metavar='LEVEL', 
        dest='verbosity',
        default='info',
        help='logging verbosity level [default: %(default)s]', choices=log_levels.keys()
    )
    
    # application version
    p.add_argument('--version', action='version', version='%(prog)s 0.1')
    
    # -- sub parsers for each action --
    s = p.add_subparsers(dest='action')
    c = s.add_parser( 'populate', help='populate samples for library',
        description='match each read in file to regions with igblast and store results in the library. takes read from stdin'
    )
    c.add_argument(
        '-l', '--library',
        dest='library', 
        metavar='NAME', 
        required=True,
        help='library name to to put results in')
        
    c = s.add_parser( 'view', help='view alignment for samples',
        description='display an alignment for samples'
    )
    c.add_argument(
        '-l', '--library',
        dest='library', 
        metavar='NAME', 
        required=True,
        help='constraint library name')
        
    c.add_argument(
        '-d', '--id',
        dest='id', 
        metavar='ID', 
        help='constraint sample id')
        
    c = s.add_parser( 'info', help='view JSON records for samples',
        description='display the JSON record of the each sample'
    )
    c.add_argument(
        '-l', '--library',
        dest='library', 
        metavar='NAME', 
        required=True,
        help='constraint library name')
        
    c.add_argument(
        '-d', '--id',
        dest='id', 
        metavar='ID', 
        help='constraint sample id')
        
    c = s.add_parser( 'simulate', help='simulate sample analysis',
        description='analyze records and print alignment and JSON records without saving to the database'
    )
    c.add_argument(
        '-l', '--library',
        dest='library', 
        metavar='NAME', 
        required=True,
        help='constraint library name')

    c.add_argument(
        '-j', '--json',
        dest='json', 
        action='store_true',
        help='emit json info')

    c.add_argument(
        '-a', '--alignment',
        dest='alignment', 
        action='store_true',
        help='emit alignment diagram')

    c = s.add_parser( 'load-reference', help='load imgt reference fasta',
        description='Load reference sequences from imgt fasta files'
    )
    c.add_argument('path', metavar='PATH', nargs='*', help='list of files to load')
    
    c = s.add_parser( 'to-blast-fasta', help='dump a fasta file for igblast',
        description='dump all samples for the region into a fasta file formatted for igblast'
    )
    c.add_argument(
        '-r', '--region',
        dest='region', 
        metavar='REGION', 
        choices=expression['regions'],
        required=True,
        help='region to dump [ {} ]'.format(', '.join(expression['regions'])))
        
    c = s.add_parser( 'to-auxiliary', help='dump auxiliary file for igblast',
        description='dump all samples for the region into an igblast auxiliary file'
    )
    c.add_argument(
        '-r', '--region',
        dest='region', 
        metavar='REGION', 
        choices=expression['regions'],
        required=True,
        help='region to dump [ {} ]'.format(', '.join(expression['regions'])))
        
    c = s.add_parser( 'rebuild', help='rebuild indexes',
        description='rebuild database indexes'
    )
    
    for k,v in vars(p.parse_args()).items():
        if v is not None:
            env[k] = v
            
    if 'action' not in env:
        p.print_help()
        
    return env

def decode_query(env):
    query = {}
    for k in [
        'library',
        'id',
        'gapped',
        'valid',
        'in frame',
        'premature termination'
    ]:
        if k in env:
            query[k] = env[k]
    return query

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    
    env = decode_cli()
    logging.getLogger().setLevel(log_levels[env['verbosity']])
    
    if 'action' in env:
        pipeline = Pipeline()
        if env['action'] == 'load-reference':
            for path in env['path']:
                pipeline.populate_reference(path)
                
        elif env['action'] == 'to-blast-fasta':
            pipeline.reference_to_blast_fasta(env['region'])
            
        elif env['action'] == 'to-auxiliary':
            pipeline.reference_to_auxiliary(env['region'])
                
        elif env['action'] == 'populate':
            pipeline.populate(env['library'])
            
        elif env['action'] == 'view':
            pipeline.view(decode_query(env))
            
        elif env['action'] == 'info':
            pipeline.info(decode_query(env))
            
        elif env['action'] == 'simulate':
            pipeline.simulate(env['library'], env['json'], env['alignment'])
            
        elif env['action'] == 'rebuild':
            pipeline.rebuild()
            
        pipeline.close()


if __name__ == '__main__':
    main()


