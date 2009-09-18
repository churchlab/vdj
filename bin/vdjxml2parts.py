#! /usr/bin/env python

import sys
import optparse

import vdj

parser = optparse.OptionParser()
parser.add_option('-b','--basename')
(options, args) = parser.parse_args()

if len(args) == 1:
    inhandle = open(args[0],'r')
    outhandle = sys.stdout
else:
    raise Exception, "Must provide at least a base outputname"

vdj.vdjxml2fasta(inhandle,outhandle)



if len(args) == 1:
    inhandle = open(args[0],'r')
elif len(args) == 0:
    inhandle = sys.stdin
else:
    raise Exception, "Too many arguments."

cmd1 = 'fasta2vdjxml.py'
cmd2 = 'size_select.py --min %d --max %d' % (options.min,options.max)
cmd3 = 'barcode_id.py --barcodes %s' % (options.barcodes_fasta)
cmd = ' | '.join([cmd1,cmd2,cmd3])