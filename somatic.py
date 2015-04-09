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
        self._reverse = None
        if node is None:
            self.node = {
                'nucleotide': None, 
                'quality': None,
                'read frame': 0,
                'strand': True
            }
        else:
            self.node = node

    def reset(self):
        for k in ['codon']:
            if k in self.node:
                del self.node[k]
        self._reverse = None

    def clone(self):
        return Sequence({
            'nucleotide': self.nucleotide, 
            'quality': self.quality,
            'read frame': self.read_frame,
            'strand': self.strand
        })

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
        return self.node['nucleotide']

    @nucleotide.setter
    def nucleotide(self, value):
        self.node['nucleotide'] = value
        self.reset()

    @property
    def quality(self):
        if 'quality' not in self.node:
            self.node['quality'] = None
        return self.node['quality']

    @quality.setter
    def quality(self, value):
        self.node['quality'] = value
        self.reset()

    @property
    def codon(self):
        if 'codon' not in self.node:
            if self.nucleotide:
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
                self.node['codon'] = ''.join(codon)
            else:
                self.node['codon'] = ''
        return self.node['codon']

    @property
    def length(self):
        return len(self.nucleotide)

    @property
    def reverse(self):
        if self._reverse is None:
            self._reverse = Sequence({
                'nucleotide': ''.join([ expression['complement'][b] for b in list(self.nucleotide) ][::-1]), 
                'quality': ''.join(self.quality[::-1]),
                'read frame': (self.length - self.read_frame) % 3,
                'strand': not self.strand
            })
        return self._reverse


class Reference(object):
    def __init__(self):
        self.log = logging.getLogger('Reference')
        self._connection = None

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

    def populate(self, path):
        sample = None
        collection = self.database['reference']
        with io.open(path, 'rb') as fasta:
            for line in fasta:
                if line:
                    line = line.strip().decode('utf8')
                    if line[0] == '>':
                        # at the start of a new record save the completed one
                        if sample is not None and sample['sequence']:
                            collection.save(sample)
                            
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
            if sample is not None and sample['sequence']:
                collection.save(sample)

    def to_blast_fasta(self, region):
        collection = self.database['reference']
        query = { 'region': region, 'strand': True }
        cursor = collection.find(query)
        for sample in cursor:
            print('>{}'.format(sample['allele name']))
            length = len(sample['sequence'])
            begin = 0
            end = 0
            while end < length:
                end = min(begin + expression['fasta line length'], length)
                print(sample['sequence'][begin:end])
                begin = end

    def to_auxiliary(self, region):
        collection = self.database['reference']
        query = { 'region': region }
        cursor = collection.find(query)
        for sample in cursor:
            if 'read frame' in sample:
                print('{}\t{}\t{}H'.format(
                    sample['allele name'], 
                    sample['read frame'], 
                    sample['region'])
                )


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
    def valid(self):
        return self.node['valid']

    @valid.setter
    def valid(self, value):
        self.node['valid'] = value

    def reverse(self):
        if self.sequence:
            self.sequence = self.sequence.reverse
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
            hit['framed'] = reference['framed']
            reference_frame = reference['read frame'] if hit['framed'] else 0
            
            if 'strain' in reference:
                hit['strain'] = reference['strain']
                
            if 'functionality' in reference:
                if 'F' in reference['functionality']:
                    hit['functionality'] = 'F'
                if 'P' in reference['functionality']:
                    hit['functionality'] = 'P'
                if 'ORF' in reference['functionality']:
                    hit['functionality'] = 'O'
                
            hit['subject'] = Sequence({
                'nucleotide': reference['sequence'][hit['subject start']:hit['subject end']],
                'strand': reference['strand'],
                'read frame': 2 - (hit['subject start'] - reference_frame - 1) % 3
            })
            
            hit['query'] = Sequence({
                'nucleotide': self.sequence.nucleotide[hit['query start']:hit['query end']],
                'strand': self.sequence.strand,
                'read frame': hit['subject'].read_frame
            })
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
                hit['overlap'] = 0
                hit['exclusive query start'] = hit['query start']
                hit['exclusive query end'] = hit['query end']
                if self.region['VH']['query end'] > hit['query start']:
                    hit['exclusive query start'] = self.region['VH']['query end']
                    hit['overlap'] += (hit['exclusive query start'] - hit['query start'])
                    
                if hit['query end'] > self.region['JH']['query start']:
                    hit['exclusive query end'] = self.region['JH']['query start']
                    hit['overlap'] += (hit['query end'] - hit['exclusive query end'])
                    
                hit['exclusive score'] = hit['bit score'] - hit['overlap'] * expression['d overlap penalty factor']
                hit['exclusive alignment length'] = hit['exclusive query end'] - hit['exclusive query start']
                if hit['exclusive query end'] > hit['exclusive query start']:
                    start = hit['exclusive query start'] - hit['query start']
                    end = start + hit['exclusive alignment length']
                    similar = 0
                    different = 0
                    longest = 0
                    stretch = 0
                    for index in range(start, end):
                        if hit['subject'].nucleotide[index] == hit['query'].nucleotide[index]:
                            similar += 1
                            stretch += 1
                            longest = max(longest, stretch)
                        else:
                            different += 1
                            stretch = 0
                    hit['exclusive identical'] = float(similar) / float(hit['exclusive alignment length'])
                    if (longest >= expression['minimum d alignment'] and
                        hit['exclusive identical'] >= expression['minimum d identity']):
                        possible.append(hit)
        if possible:
            possible.sort(reverse=True, key=lambda x: x['exclusive score'])
            hit = possible[0]
            self.region['DH'] = hit
            hit['picked'] = True
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
            

    def analyze(self):
        if not self.gapped:
            valid = self._complement_from_reference()
            valid = valid and self._pick_jh_region()
            valid = valid and self._assign_frame()
            valid = valid and self._pick_vh_region()
            valid and self._pick_dh_region()
            valid and self._identify_cdr3()
            self.valid = valid
        else:
            self.log.info('Skipping sample with gapped alignment %s', self.id)
            self.valid = False

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
                    
            if self.region['VH'] and self.region['JH']:
                if self.region['VH']['in frame'] == 'Y' and self.region['JH']['in frame'] == 'Y':
                    buffer.write('in frame')
                else:
                    buffer.write('out of frame')
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
            if not sample.gapped:
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

    def fill(self):
        if self.read():
            buffer = self.search()
            if buffer:
                self.parse_igblast(buffer)
                for sample in self.buffer:
                    # sample.compress()
                    sample.analyze()
        return not self.empty

    def read(self):
        self.node = {
            'buffer': [],
            'document': [],
            'lookup': {}
        }
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


class Pipeline(object):
    def __init__(self):
        self.log = logging.getLogger('Pipeline')
        self.count = 0
        self._connection = None

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

    def reference_for(self, allele_name):
        return self.database['reference'].find_one({'allele name': allele_name})

    def rebuild(self):
        index = [
            {
                'collection': 'reference',
                'index': [
                    { 'key': [('region', ASCENDING )], 'unique': False, 'name': 'region' },
                    { 'key': [( 'allele name', ASCENDING )], 'unique': True, 'name': 'allele name' },
                    { 'key': [( 'subgroup', ASCENDING )], 'unique': False, 'name': 'subgroup' },
                    { 'key': [( 'gene name', ASCENDING )], 'unique': False, 'name': 'gene name' }
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

    def populate(self, library):
        block = Block(self)
        collection = self.database[library]
        while block.fill():
            print(str(block))
            print(to_json(block.document))
            # result = collection.insert_many(block.document)
            # self.log.debug('%s so far', self.count)
            self.count += block.size
            
    def view(self, library=None, id=None):
        block = Block(self)
        while block.fill():
            print(str(block))
            print(to_json(block.document))
            self.count += block.size


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
        help='constraint library name')
        
    c.add_argument(
        '-d', '--id',
        dest='id', 
        metavar='ID', 
        help='constraint sample id')
        
    c = s.add_parser( 'rebuild', help='rebuild indexes',
        description='rebuild database indexes'
    )
    
    for k,v in vars(p.parse_args()).items():
        if v is not None:
            env[k] = v
            
    if 'action' not in env:
        p.print_help()
        
    return env

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    
    env = decode_cli()
    logging.getLogger().setLevel(log_levels[env['verbosity']])
    
    if 'action' in env:
        if env['action'] == 'load-reference':
            reference = Reference()
            for path in env['path']:
                reference.populate(path)
                
        elif env['action'] == 'to-blast-fasta':
            reference = Reference()
            reference.to_blast_fasta(env['region'])
            
        elif env['action'] == 'to-auxiliary':
            reference = Reference()
            reference.to_auxiliary(env['region'])
                
        elif env['action'] == 'populate':
            pipeline = Pipeline()
            pipeline.populate(env['library'])
            
        elif env['action'] == 'view':
            pipeline = Pipeline()
            pipeline.view()
            # pipeline.view(env['library'], env['id'])
            
        elif env['action'] == 'rebuild':
            pipeline = Pipeline()
            pipeline.rebuild()

if __name__ == '__main__':
    main()


