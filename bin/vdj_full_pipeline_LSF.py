#! /usr/bin/env python

import os
import sys
import optparse
import subprocess

from Bio import SeqIO
from Bio.Alphabet import generic_dna

import lsf
import seqtools
import vdj
import vdj.pipeline



# 0. SETUP
join = os.path.join
def log(s): sys.stdout.write(s); sys.stdout.flush()

option_parser = optparse.OptionParser()
(options,args) = option_parser.parse_args()

# Read and process jobfile
if len(args) != 1:
    raise ValueError, "Require a single input jobfile"
jobfile = args[0]
params = vdj.pipeline.parse_jobfile(jobfile)
work_dir = os.path.abspath(params['work_dir'])
basename = params['basename']

short_queue = params['short_queue']
long_queue = params['long_queue']

# Create working directories
assert not os.path.exists(work_dir)
os.mkdir(work_dir,0755)
os.mkdir(join(work_dir,'parts'),0755)
os.mkdir(join(work_dir,'logs'),0755)
os.mkdir(join(work_dir,'partitions'),0755)



# 1. SIZE SELECTION
log("Performing size selection on reads...")
min_size = int(params['min_size'])
max_size = int(params['max_size'])
size_selected_file = join(work_dir,basename + '.size%i-%i' % (min_size,max_size) + '.imgt')
with open(size_selected_file,'w') as outhandle:
    for seq in SeqIO.parse(params['input_fasta'],'fasta',generic_dna):
        if len(seq) >= min_size and len(seq) <= max_size:
            chain = vdj.ImmuneChain(seq)
            print >>outhandle, chain
log("finished\n")



# 2. SPLIT INTO PARTS
log("Splitting input into small parts...")
parts = vdj.pipeline.iterator2parts( vdj.parse_imgt(size_selected_file),
                                     join(work_dir,'parts/size_selected.imgt'),
                                     int(params['packet_size']))
log("finished\n")



# 3-7. BARCODE ID, CODING STRAND, ISOTYPE ID, VDJ CLASSIFICATION, TRANSLATION via LSF
log("Setting up LSF command...\n")
locus_options = ' '.join([' --locus %s' % locus for locus in params['locus']])
cmd = 'barcode_id.py --barcodes %s ' % params['barcode_fasta']      # 3. BARCODE IDENTIFICATION
cmd += ' | coding_strand.py' + locus_options                        # 4. CODING STRAND
if 'IGH' in params['locus']:                                        # 5. ISOTYPE ID (heavy chain only)
    cmd += ' | isotype_id.py --IGHC %s' % params['isotype_fasta']
cmd += ' | align_vdj.py' + locus_options                            # 6. VDJ CLASSIFICATION
cmd += ' | translate_chains.py'                                     # 7. TRANSLATION

# submit cmd to LSF for each part
log("Submitting jobs to LSF...")
jobIDs = []
logfiles = []
outnames = []
for part in parts:
    partID = part.split('.')[-1]
    partoutname = join(work_dir,'parts/aligned.imgt.'+partID)
    outnames.append(partoutname)
    curr_cmd = 'cat %s | ' + cmd + ' > %s'
    curr_cmd = curr_cmd % (part,partoutname)
    logfile = join(work_dir,'logs/alignment.log.'+partID)
    jobID = lsf.submit_to_LSF(short_queue,logfile,curr_cmd)
    logfiles.append(logfile)
    jobIDs.append(jobID)
log("finished\n")

log("Waiting for LSF jobs...")
lsf.wait_for_LSF_jobs(jobIDs,logfiles)
log("finished\n")



# 8. CONCAT PARTS
log("Concatenating pieces...")
aligned_file = join(work_dir,basename+'.aligned.imgt')
cmd = 'find %s -name "aligned.imgt.*" | xargs cat > %s' % (join(work_dir,'parts'),aligned_file)
p = subprocess.Popen(cmd,shell=True)
p.wait()
log("finished\n")



# 9. FILTER VJ
log("Filtering for chains with VJ alignments...")
vj_filtered_file = join(work_dir,basename+'.VJ_filtered.imgt')
with open(vj_filtered_file,'w') as outhandle:
    for chain in vdj.parse_imgt(aligned_file):
        if hasattr(chain,'v') and hasattr(chain,'j'):
            print >>outhandle, chain
log("finished\n")



# ###########################################################################

# =======
# = OLD =
# =======

# # 10. PARTITION VJ
# sys.stdout.write("Partitioning data by VJ alignment..."); sys.stdout.flush()
# inhandle = open(join(analysis_dir,vj_filtered_file),'r')
# partition_files = vdj.pipeline.partition_VJ(inhandle,join(work_dir,partition_dir,basename))
# inhandle.close()
# sys.stdout.write("finished\n"); sys.stdout.flush()
# 
# # 11. CLUSTER CDR3s
# sys.stdout.write("Setting up clustering jobs...\n"); sys.stdout.flush()
# cmd = 'cluster_cdr3.py --cutoff %f --linkage %s' % (clustering_cutoff,clustering_linkage)
# 
# # submit LSF
# sys.stdout.write("Submitting clustering jobs to LSF..."); sys.stdout.flush()
# jobIDs = []
# logfiles = []
# outnames = []
# for partition in partition_files:
#     partitionoutname = '.'.join(partition.split('.')[:-1])+'.clustered.vdjxml'
#     outnames.append(partitionoutname)
#     tag = os.path.splitext(os.path.basename(partition))[0].split('.')[-1]
#     curr_cmd = 'cat %s | ' + cmd + ' --tag %s' + ' > %s'
#     curr_cmd = curr_cmd % (partition,tag,partitionoutname)
#     logfile = join(work_dir,log_dir,'clustering.')+tag+'.log'
#     jobID = lsf.submit_to_LSF(long_queue,logfile,curr_cmd)
#     logfiles.append(logfile)
#     jobIDs.append(jobID)
# sys.stdout.write("finished\n"); sys.stdout.flush()
# 
# sys.stdout.write("Waiting for LSF jobs..."); sys.stdout.flush()
# lsf.wait_for_LSF_jobs(jobIDs,logfiles)
# sys.stdout.write("finished\n"); sys.stdout.flush()
# 
# # 12. CONCAT PARTITIONS
# sys.stdout.write("Concatenating partitions..."); sys.stdout.flush()
# outhandle = open(join(analysis_dir,clustered_file),'w')
# vdj.pipeline.cat_vdjxml(outnames,outhandle)
# outhandle.close()
# sys.stdout.write("finished\n"); sys.stdout.flush()
