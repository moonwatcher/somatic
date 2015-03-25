#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import logging
import re
import json
import io

from subprocess import Popen, PIPE
from io import StringIO
from pymongo import MongoClient
from argparse import ArgumentParser

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
    'regions': ['V', 'D', 'J'],
    'fasta line length': 80,
    'igblast result start': re.compile('^# IGBLASTN '),
    'igblast query': re.compile('^# Query: (P:<query_id>.*)$'),
    'igblast hit table': re.compile('^# Hit table'),
    'igblast result end': re.compile('^# BLAST processed '),
    'nucleotide sequence': re.compile('^[atcgnATCGN]+$'),
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
    def __init__(self, node=None):
        self.log = logging.getLogger('Sample')
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


class Processor(object):
    def __init__(self):
        self.log = logging.getLogger('Processor')
        self.buffer_size = 128
        self.buffer = []
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

    @property
    def full(self):
        return not len(self.buffer) < self.buffer_size

    def populate(self, library):
        collection = self.database[library]
        while self.read_block():
            self.search_igblast()
            self.compress()
            
            batch = [ s.node for s in self.buffer ]
            print(to_json(batch))
            # result = collection.insert_many(batch)
            self.count += len(self.buffer)
            self.log.debug('%s so far', self.count)

    def compress(self):
        for sample in self.buffer:
            m = []
            for match in sample.match:
                try:
                    m.append({ 'hit': expression['igblast compressed hit'].format(**match) })
                except KeyError as e:
                    self.log.error(e)
                    self.log.error(match)
            sample.match.clear()
            sample.match.extend(m)

    def read_block(self):
        state = None
        self.buffer = []
        self.lookup = {}
        sample = Sample()
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
                    self.buffer.append(sample)
                    self.lookup[sample.id] = sample
                    state = None
                    
                    if self.full:
                        break
                    else:
                        sample = Sample()
            else:
                break
                
        return len(self.buffer) > 0

    def to_fasta(self):
        fasta = []
        for sample in self.buffer:
            fasta.append('>{}'.format(sample.id))
            begin = 0
            end = 0
            while end < sample.length:
                end = min(begin + expression['fasta line length'], sample.length)
                fasta.append(sample.sequence[begin:end])
                begin = end
        return fasta

    def search_igblast(self):
        fasta = self.to_fasta()
        process = Popen(
            args=expression['igblast command'],
            cwd='/Users/lg/code/somatic/igblast',
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate(input='\n'.join(fasta).encode('utf8'))
        if output:
            buffer = StringIO(output.decode('utf8'))
            self.parse_igblast_output(buffer)

    def parse_igblast_output(self, buffer):
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
        processor = Processor()
        processor.populate(env['library'])

if __name__ == '__main__':
    main()

    