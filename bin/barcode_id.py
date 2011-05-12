#! /usr/bin/env python

import sys
import optparse

import seqtools

import vdj
import vdj.pipeline

parser = optparse.OptionParser()
parser.add_option('-b','--barcodes',dest='barcodes_fasta')
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

# NOTE: all barcodes must be the same length

barcodes = vdj.pipeline.load_barcodes(options.barcodes_fasta)

# iterate through chains
for chain in vdj.parse_imgt(inhandle):
    vdj.pipeline.id_barcode(chain,barcodes)
    print >>outhandle, chain
