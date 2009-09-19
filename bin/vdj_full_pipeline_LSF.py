#! /usr/bin/env python

import subprocess
import sys
import optparse

import vdj

parser = optparse.OptionParser()
parser.add_option('-q','--queue')
parser.add_option('-o','--LSFoutput')
parser.add_option('-p','--packetsize',type='int')
parser.add_option('-b','--barcodes',dest='barcodes_fasta')
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
cmd3 = 'vdjxml2parts.py --packetsize %d --basename %s' % (options.packetsize,args[0])
cmd = ' | '.join([cmd1,cmd2,cmd3])
p = subprocess.Popen(cmd,shell=True,stdin=inhandle,stdout=subprocess.PIPE)
parts = [f.strip() for f in p.stdout.readlines()]
outparts = [part+'.out' for part in parts]

cmd4 = 'cat %s'
cmd5 = 'barcode_id.py --barcodes %s' % (options.barcodes_fasta)
cmd6 = 'positive_strand.py'
cmd7 = 'align_vdj.py > %s'
cmd = ' | '.join([cmd4,cmd5,cmd6,cmd7])
jobs = []
for part,outpart in zip(parts,outparts):
    jobID = vdj.LSF.submit_to_LSF(options.queue,options.LSFoutput,cmd % (part,outpart))
    jobs.append(jobID)
vdj.LSF.wait_for_LSF_jobs(jobs)

file_list = ' '.join(outparts)
subprocess.Popen('cat ' + file_list,shell=True,stdout=outhandle)
