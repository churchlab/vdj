#! /usr/bin/env python

import subprocess
import sys
import optparse

import vdj

parser = optparse.OptionParser()
parser.add_option('-b','--barcodes',dest='barcodes_fasta')
parser.add_option('-c','--cutoff',default=4.5,type='float')
parser.add_option('-i','--IGHC',dest='ighc_fasta')
parser.add_option('-m','--min',type='int')
parser.add_option('-M','--max',type='int')
(options, args) = parser.parse_args()

if len(args) == 2:
    inhandle = open(args[0],'r')
    outhandle = open(args[1],'w')
else:
    raise Exception, "Must have an input file and output file."

cmd1 = 'fasta2vdjxml.py'
cmd2 = 'size_select.py --min %d --max %d' % (options.min,options.max)
cmd3 = 'barcode_id.py --barcodes %s' % (options.barcodes_fasta)
cmd = ' | '.join([cmd1,cmd2,cmd3])
p = subprocess.Popen(cmd,shell=True,stdin=inhandle,stdout=subprocess.PIPE)
