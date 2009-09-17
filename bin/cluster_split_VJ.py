#! /usr/bin/env python

import sys
import optparse
import subprocess

import vdj

parser = optparse.OptionParser()
parser.add_option('-c','--cutoff',default=4.5,type='float')
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

(VJ_parts,VJ_IDs) = vdj.split_vdjxml_into_VJ_parts(inhandle,outname)

VJ_parts_clustered = []
for (vj_file,vj_id) in zip(VJ_parts,VJ_IDs):
    vj_file_clustered = vj_file + '.clustered'
    VJ_parts_clustered.append(vj_file_clustered)
    params = {'cutoff':options.cutoff,
                 'tag':options.tag,
              'infile':vj_file,
             'outfile':vj_file_clustered}
    cluster_cmd = r'python cluster_cdr3 --cutoff %(cutoff)f --tag %(tag)s %(infile)s %(outfile)s' % params
    p = subprocess.Popen(cluster_cmd,shell=True)
    # perform serially:
    # vdj.cluster_chains(options.cutoff,vj_id,vj_file,vj_file_clustered)

# WHEN ALL JOBS ARE DONE

for chain in vdj.parse_VDJXML_parts(VJ_parts_clustered):
    print >>outhandle, chain