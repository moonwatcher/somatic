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
from io import StringIO

def usage():
    print('usage: path int:start int:end')

def crop(path, start, end):
    reference = StringIO()
    with io.open(path, 'r') as file:
        try:
            for line in file:
                if end is not None and reference.tell() > end:
                    break
                    
                line = line.strip()
                if line[0] != '>':
                    reference.write(line)
                    
            reference.seek(start)
            if end is not None:
                print(reference.read(end - start))
            else:
                print(reference.read())
                
        except OSError as e:
            print('{} {}'.format(e.strerror, path))

if len(sys.argv) > 2:
    path = sys.argv[1]
    try:
        start = int(sys.argv[2])
        if len(sys.argv) > 3:
            end = int(sys.argv[3])
        else:
            end = None
    except ValueError as e:
        usage()
    else:
        crop(path, start, end)
else:
    usage()
