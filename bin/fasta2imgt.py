#! /usr/bin/env python

import sys
import optparse

from Bio import SeqIO

import vdj

parser = optparse.OptionParser()
(options, args) = parser.parse_args()

if len(args) == 2:
    inhandle = open(args[0],'r')
    outhandle = open(args[1],'w')
elif len(args) == 1:
    inhandle = open(args[0],'r')
    outhandle = sys.stdout
elif len(args) == 0:
    inhandle = sys.stdin
    outhandle = sys.stdout
else:
    raise Exception, "Wrong number of arguments."

for record in SeqIO.parse(inhandle,'fasta'):
    chain = vdj.ImmuneChain(record)
    print >>outhandle, chain
