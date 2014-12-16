#! /usr/bin/env python
# Copyright 2014 Uri Laserson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import argparse

from Bio import SeqIO
from Bio.Alphabet import generic_dna

import vdj

argparser = argparse.ArgumentParser(description=None)
argparser.add_argument('input',nargs='?',type=argparse.FileType('r'),default=sys.stdin)
argparser.add_argument('output',nargs='?',type=argparse.FileType('w'),default=sys.stdout)
args = argparser.parse_args()

for record in SeqIO.parse(args.input,'fasta',generic_dna):
    chain = vdj.ImmuneChain(record.upper())
    print >>args.output, chain
