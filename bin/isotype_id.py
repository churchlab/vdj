#! /usr/bin/env python

import sys
import warnings
import optparse

import seqtools

import vdj

parser = optparse.OptionParser()
parser.add_option('-i','--IGHC',dest='ighc_fasta')
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

# load isotypes
ighcip = open(options.ighc_fasta,'r')
isotypes = {}
for (descr,seq) in seqtools.FastaIterator(ighcip):
    isotypes[seq.upper()] = record.id
ighcip.close()

for chain in parse_VDJXML(inhandle):
    if not chain.has_tag('positive') and not chain.has_tag('coding'):
        warnings.warn('chain %s may not be the correct strand' % chain.descr)
    
    for iso in isotypes.iteritems():
        if iso[0] in chain.seq[-50:]:   # arbitrary cutoff from 3' end
            chain.c = iso[1]
    print >>outhandle, chain
