#! /usr/bin/env python

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
