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
    return json.dumps(node, sort_keys=False, indent=4)

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
    'regions': ['V', 'D', 'J'],
    'fasta line length': 80,
    'igblast result start': re.compile('^# IGBLASTN '),
    'igblast query': re.compile('^# Query: (P:<query_id>.*)$'),
    'igblast hit table': re.compile('^# Hit table'),
    'igblast result end': re.compile('^# BLAST processed '),
    'nucleotide sequence': re.compile('^[atcgnATCGN]+$'),
    'expand igblast hit': re.compile(
        r"""
        (?P<segment>[VDJ]),
        (?P<subject_id>[^,]+),
        (?P<query_start>[^,]+),
        (?P<query_end>[^,]+),
        (?P<subject_start>[^,]+),
        (?P<subject_end>[^,]+),
        (?P<gap_openings>[^,]+),
        (?P<gaps>[^,]+),
        (?P<mismatches>[^,]+),
        (?P<percentage_identical>[^,]+),
        (?P<bit_score>[^,]+),
        (?P<evalue>[^,]+),
        (?P<alignment_length>[^,]+),
        (?P<subject_strand>[^,]),
        (?P<query_strand>[^,])
        """,
        re.VERBOSE
    ),
    'igblast hit': re.compile(
        r"""
        (?P<segment>[VDJ])\t
        (?:(?P<query_polarity>reversed)\|)?(?:[^,]+)\t
        (?P<subject_id>[^,]+)\t
        (?P<query_start>[^,]+)\t
        (?P<query_end>[^,]+)\t
        (?P<subject_start>[^,]+)\t
        (?P<subject_end>[^,]+)\t
        (?P<gap_openings>[^,]+)\t
        (?P<gaps>[^,]+)\t
        (?P<mismatches>[^,]+)\t
        (?P<percentage_identical>[^,]+)\t
        (?P<bit_score>[^,]+)\t
        (?P<evalue>[^,]+)\t
        (?P<alignment_length>[^,]+)\t
        (?P<subject_strand>[^,]+)
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
    'igblast compressed hit':'{segment},{subject_id},{query_start},{query_end},{subject_start},{subject_end},{gap_openings},{gaps},{mismatches},{percentage_identical},{bit_score},{evalue},{alignment_length},{subject_strand},{query_strand}',
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
        (?P<reading_frame>[123]|NR)\|
        [^|]*\|
        [^|]*\|
        [^|]*\|
        [^|]*\|
        [^|]*\|
        [^|]*\|
        (?:rev-compl)?\s*\|$
        """,
        re.VERBOSE
    ),
}

    
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
                            for k,v in match.groupdict().items():
                                if v:
                                    try:
                                        if k in ['start', 'end', 'length']:
                                            v = int(v)
                                        elif k == 'reading_frame':
                                            if v == 'NR':
                                                v = None
                                            else:
                                                v = int(v)
                                    except ValueError as e:
                                        self.log.error('unable to decode %s as int for %s', v, k)
                                    else:
                                        if v is not None:
                                            sample[k] = v
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
        query = { 'region': region }
        cursor = collection.find(query)
        for sample in cursor:
            print('>{}'.format(sample['allele_name']))
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
            if 'reading_frame' in sample:
                print('{}\t{}\t{}H'.format(
                    sample['allele_name'], 
                    sample['reading_frame'] - 1, 
                    sample['region'])
                )


class Sample(object):
    def __init__(self, pipeline, node=None):
        self.log = logging.getLogger('Sample')
        self.pipeline = pipeline
        self.node = node
        
        if self.node is None:
            self.node = {}

    @property
    def id(self):
        return self.node['id']

    @id.setter
    def id(self, value):
        self.node['id'] = value

    @property
    def sequence(self):
        return self.node['sequence']

    @sequence.setter
    def sequence(self, value):
        self.node['sequence'] = value

    @property
    def quality(self):
        return self.node['sequence']

    @quality.setter
    def quality(self, value):
        self.node['quality'] = value

    @property
    def length(self):
        return len(self.node['sequence'])

    @property
    def match(self):
        if 'match' not in self.node:
            self.node['match'] = []
        return self.node['match']

    @property
    def codon_sequence(self):
        if 'codon_sequence' not in self.node:
            if self.reading_frame and self.sequence:
                self.node['codon_sequence'] = self.to_codon(self.sequence, self.reading_frame)
            return self.node['codon_sequence']
            
        if 'codon_sequence' not in self.node:
            self.node['codon_sequence'] = None
            
        return self.node['codon_sequence']

    @property
    def reverse_complement_codon_sequence(self):
        if 'reverse_complement_codon_sequence' not in self.node:
            if self.reading_frame and self.reverse_complement_sequence:
                self.node['reverse_complement_codon_sequence'] = self.to_codon(self.reverse_complement_sequence, self.reading_frame)
            return self.node['reverse_complement_codon_sequence']
            
        if 'reverse_complement_codon_sequence' not in self.node:
            self.node['reverse_complement_codon_sequence'] = None
            
        return self.node['reverse_complement_codon_sequence']

    @property
    def reverse_complement_sequence(self):
        if 'reverse complement sequence' not in self.node:
            self.node['reverse complement sequence'] = ''.join([ expression['complement'][b] for b in list(self.sequence) ][::-1])
        return self.node['reverse complement sequence']

    @property
    def reverse_complement_quality(self):
        if 'reverse complement quality' not in self.node:
            self.node['reverse complement quality'] = ''.join(self.quality[::-1])
        return self.node['reverse complement sequence']

    @property
    def reading_frame(self):
        if 'reading_frame' not in self.node:
            if 'match' in self.node and self.node['match']:
                if 'reading_frame' in self.node['match'][0]:
                    self.node['reading_frame'] = self.node['match'][0]['reading_frame']
                    
        if 'reading_frame' not in self.node:
            self.node['reading_frame'] = None
        
        return self.node['reading_frame']

    def to_codon(self, sequence, offset):
        length = len(sequence)
        start = offset - 1
        end = start + 3
        amino = []
        while(end <= length):
            acid = sequence[start:end]
            if 'N' in acid:
                amino.append('X')
            else:
                amino.append(expression['nucleic to amino'][acid])
            start = end
            end = start + 3
        return ''.join(amino)

    def expand(self):
        if 'match' in self.node:
            hits = []
            for r in self.node['match']:
                if 'hit' in r:
                    match = expression['expand igblast hit'].search(r['hit'])
                    if match:
                        hit = match.groupdict()
                        hit['hit'] = r['hit']
                        for k in [
                            'query_start',
                            'query_end',
                            'subject_start',
                            'subject_end',
                            'gap_openings',
                            'gaps',
                            'mismatches',
                            'alignment_length'
                        ]:
                            try: hit[k] = int(hit[k])
                            except ValueError as e:
                                self.log.warning('could not parse value %s for %s as int', hit[k], k)
                                hit[k] = None
                                
                        for k in [
                            'percentage_identical',
                            'bit_score',
                            'evalue',
                        ]: 
                            try: hit[k] = float(hit[k])
                            except ValueError:
                                self.log.warning('could not parse value %s for %s as int', hit[k], k)
                                hit[k] = None
                        hits.append(hit)
        self.node['match'] = hits

    def compress(self):
        for match in self.match:
            try:
                match['hit'] = expression['igblast compressed hit'].format(**match)
            except KeyError as e:
                self.log.error(e)
                self.log.error(match)

    def analyze(self):
        for match in self.match:
            reference = self.pipeline.reference_for(match['subject_id'])
            if 'reading_frame' in reference:
                match['subject_sequence'] = reference['sequence'][match['subject_start'] - 1:match['subject_end']]
                if match['subject_strand'] == match['query_strand']:
                    match['query_sequence'] = self.sequence[match['query_start'] - 1:match['query_end']]
                else:
                    match['query_sequence'] = self.reverse_complement_sequence[match['query_start'] - 1:match['query_end']]
                match['reading_frame'] = 3 - (match['subject_start'] - reference['reading_frame'] - 1) % 3
                match['amino_query_sequence'] = self.to_codon(match['query_sequence'], match['reading_frame'])
                match['amino_subject_sequence'] = self.to_codon(match['subject_sequence'], match['reading_frame'])

    def to_json(self):
        def handler(o):
            result = None
            if isinstance(o, datetime):
                result = o.isoformat()
            if isinstance(o, ObjectId):
                result = str(o)
            if isinstance(o, set):
                result = list(o)
            return result

        return json.dumps(self.node, sort_keys=False, ensure_ascii=False, indent=4, default=handler)

    def __str__(self):
        buffer = StringIO()
        if self.reading_frame:
            # print the coordinate system
            if self.reading_frame > 1:
                buffer.write('{: <{}} '.format(1, self.reading_frame - 1))
            for index in range(self.reading_frame, self.length, 3):
                buffer.write('{: <3} '.format(index))
            buffer.write('\n')
            
            # print the sample nucleotide sequence
            if self.reading_frame > 1:
                buffer.write(self.reverse_complement_sequence[0:self.reading_frame - 1])
                buffer.write(' ')
            for index in range(self.reading_frame - 1, self.length, 3):
                buffer.write(self.reverse_complement_sequence[index:index + 3])
                buffer.write(' ')
            buffer.write('\n')
            
            # print the sample codon sequence
            if self.reading_frame > 1:
                buffer.write('{: <{}} '.format(' ', self.reading_frame - 1))
            for codon in self.reverse_complement_codon_sequence:
                buffer.write(codon)
                buffer.write('   ')
            buffer.write('\n')
            
            for match in self.match:
                if match['segment'] == 'D':
                    reading_frame = 3 - (match['query_start'] - self.reading_frame - 1) % 3
                    amino_subject_sequence = self.to_codon(match['subject_sequence'], reading_frame)
                else:
                    reading_frame = match['reading_frame']
                    amino_subject_sequence = match['amino_subject_sequence']
                    
                # print the nucleotide sequence
                mask = []
                for index in range(match['query_start'] - 1, match['query_end']):
                    if self.reverse_complement_sequence[index] == match['subject_sequence'][index - match['query_start'] + 1]:
                        mask.append('-')
                    else:
                        mask.append(match['subject_sequence'][index - match['query_start'] + 1])
                mask = ''.join(mask)
                    
                if match['query_start'] > 1:
                    buffer.write('{: <{}} '.format(' ', int((match['query_start'] - self.reading_frame) / 3) - 1 + match['query_start']))
                    
                if reading_frame > 1:
                    buffer.write(mask[0:reading_frame - 1])
                    buffer.write(' ')
                    
                for index in range(reading_frame - 1, match['alignment_length'], 3):
                    buffer.write(mask[index:index + 3])
                    buffer.write(' ')
                buffer.write('\n')

                # print the sample codon sequence
                if match['query_start'] > 1:
                    buffer.write('{: <{}} '.format(' ', int((match['query_start'] - self.reading_frame) / 3) - 1 + match['query_start']))
                if reading_frame > 1:
                    buffer.write('{: <{}} '.format(' ', reading_frame - 1))
                
                for codon in amino_subject_sequence:
                    buffer.write(codon)
                    buffer.write('   ')
                buffer.write('\n')
                
        buffer.seek(0)
        return buffer.read()


class Block(object):
    def __init__(self, pipeline):
        self.log = logging.getLogger('Block')
        self.pipeline = pipeline
        self.reset()

    def __str__(self):
        buffer = StringIO()
        for sample in self.buffer:
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
                    sample.compress()
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
                    sample.sequence = line
                    state = 'sequence'
                    
                elif state == 'sequence':
                    if line[0] == '+':
                        state = 'redundant'
                        
                elif state == 'redundant':
                    sample.quality = line
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
            while end < sample.length:
                end = min(begin + expression['fasta line length'], sample.length)
                buffer.write(sample.sequence[begin:end])
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
            elif state == 1 and  line.startswith('# Query: '):
                id = line[9:]
                if id in self.lookup:
                    sample = self.lookup[id]
                else:
                    self.log.error('OH MY GOD!')
                state = 2
            elif state == 2 and line.startswith('# Hit table '):
                state = 3
            elif (state == 3 or state == 4) and not line.startswith('#'):
                state = 4
                match = expression['igblast hit'].search(line)
                if match:
                    hit = match.groupdict()
                    if 'query_polarity' in hit and hit['query_polarity'] == 'reversed':
                        if hit['subject_strand'] == 'plus':
                            hit['query_strand'] = 'minus'
                        else:
                            hit['query_strand'] = 'plus'
                        del hit['query_polarity']
                    else:
                        hit['query_strand'] = hit['subject_strand']
                        
                    for k in [
                        'query_strand',
                        'subject_strand'  
                    ]:
                        if hit[k] == 'plus': hit[k] = '+'
                        elif hit[k] == 'minus': hit[k] = '-'
                        else:
                            self.log.warning('could not parse value %s for %s as a strand', hit[k], k)
                            hit[k] = None
                            
                    for k in [
                        'percentage_identical',
                        'bit_score',
                        'evalue',
                    ]: 
                        try: hit[k] = float(hit[k])
                        except ValueError:
                            self.log.warning('could not parse value %s for %s as int', hit[k], k)
                            hit[k] = None
                        
                    for k in [
                        'subject_end',
                        'gap_openings',
                        'subject_start',
                        'mismatches',
                        'query_start',
                        'alignment_length',
                        'gaps',
                        'query_end',
                    ]:
                        try: hit[k] = int(hit[k])
                        except ValueError as e:
                            self.log.warning('could not parse value %s for %s as int', hit[k], k)
                            hit[k] = None
                    sample.match.append(hit)
                    
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
        return self.database['reference'].find_one({'allele_name': allele_name})
    

    def rebuild(self):
        index = [
            {
                'collection': 'reference',
                'index': [
                    { 'key': [('region', ASCENDING )], 'unique': False, 'name': 'region' },
                    { 'key': [( 'allele_name', ASCENDING )], 'unique': True, 'name': 'allele_name' },
                    { 'key': [( 'subgroup', ASCENDING )], 'unique': False, 'name': 'subgroup' },
                    { 'key': [( 'gene_name', ASCENDING )], 'unique': False, 'name': 'gene_name' }
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
            print(to_json(block.document))
            print(str(block))
            # result = collection.insert_many(block.document)
            self.count += block.size
            self.log.debug('%s so far', self.count)


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
        
    c = s.add_parser( 'rebuild', help='rebuild indexes',
        description='rebuild database indexes'
    )
    
    for k,v in vars(p.parse_args()).items():
        if v is not None:
            env[k] = v
    return env
    

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    
    env = decode_cli()
    logging.getLogger().setLevel(log_levels[env['verbosity']])
    
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
        processor = Pipeline()
        processor.populate(env['library'])

    elif env['action'] == 'analyze':
        processor = Pipeline()
        processor.analyze(env['library'])

    elif env['action'] == 'rebuild':
        processor = Pipeline()
        processor.rebuild()
        

if __name__ == '__main__':
    main()

    