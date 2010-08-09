#! /usr/bin/env python

import sys
import optparse

import seqtools

import vdj

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

# load barcodes
bcip = open(options.barcode_fasta,'r')
barcodes = {}
for (descr,seq) in seqtools.FastaIterator(bcip):
    barcodes[seq.upper()] = descr
bcip.close()

# check that barcodes meet necessary criteria
barcode_len = len(barcodes.keys()[0])
for bc in barcodes.keys():
    if len(bc) != barcode_len:
        raise Exception, "ERROR: All barcode lengths must be equal."

# iterate through chains
print >>outhandle, "<root>"
for chain in vdj.parse_VDJXML(inhandle):
    try:
        curr_barcode = barcodes[chain.seq[:barcode_len].upper()]
        chain.seq = chain.seq[barcode_len:] # prune off barcode from seq
        chain.barcode = curr_barcode
    except KeyError:    # barcode not found; print chain unchanged
        pass    # chain remains unchanged
    print >>outhandle, chain
print >>outhandle, "</root>"
