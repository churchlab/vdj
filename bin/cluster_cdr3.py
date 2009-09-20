#! /usr/bin/env python

import sys
import optparse

import vdj

parser = optparse.OptionParser()
parser.add_option('-c','--cutoff',default=4.5,type='float')
parser.add_option('-t','--tag',default='')
parser.add_option('-l','--linkage',type='choice',choices=['single','complete'],default='single')
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

vdj.cluster_chains(options.cutoff,options.tag,inhandle,outhandle,options.linkage)