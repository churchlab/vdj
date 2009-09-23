#! /usr/bin/env python

import sys
import optparse
import subprocess
import os
import tempfile

import vdj

parser = optparse.OptionParser()
parser.add_option('-c','--cutoff',default=4.5,type='float')
parser.add_option('-l','--linkage',type='choice',choices=['single','complete'],default='single')
parser.add_option('-q','--queue')
parser.add_option('-o','--LSFoutput')
(options, args) = parser.parse_args()

if len(args) == 2:
    inhandle = open(args[0],'r')
    outname = args[0]
    outhandle = open(args[1],'w')
elif len(args) == 1:
    inhandle = open(args[0],'r')
    outname = args[0]
    outhandle = sys.stdout
elif len(args) == 0:
    inhandle = sys.stdin
    outname = 'VJ_parts.vdjxml'
    outhandle = sys.stdout

print >>sys.stderr, "NOTE: chains must be filtered for valid VJ aln and junctions BEFORE clustering."

# temporary directory to dump intermediate files
tempdir = tempfile.mkdtemp(prefix=outname+'.intermediate.',dir='.')

(VJ_parts,VJ_IDs) = vdj.split_vdjxml_into_VJ_parts(inhandle,os.path.join(tempdir,outname))

VJ_parts_clustered = []
jobs = []
for (vj_file,vj_id) in zip(VJ_parts,VJ_IDs):
    vj_file_clustered = vj_file+'.clustered'    # NOTE: vj_file already has directory prefix
    VJ_parts_clustered.append(vj_file_clustered)
    params = {'cutoff':options.cutoff,
             'linkage':options.linkage,
                 'tag':vj_id,
              'infile':vj_file,
             'outfile':vj_file_clustered}
    cluster_cmd = r'cluster_cdr3.py --cutoff %(cutoff)f --tag %(tag)s --linkage %(linkage)s %(infile)s %(outfile)s' % params
    if os.stat(vj_file).st_size < 8e6:
        jobID = vdj.LSF.submit_to_LSF(options.queue,options.LSFoutput,cluster_cmd)
    else:
        jobID = vdj.LSF.submit_to_LSF(options.queue,options.LSFoutput,cluster_cmd,mem_usage=4096)
    jobs.append(jobID)

vdj.LSF.wait_for_LSF_jobs(jobs)

for chain in vdj.parse_VDJXML_parts(VJ_parts_clustered):
    print >>outhandle, chain
