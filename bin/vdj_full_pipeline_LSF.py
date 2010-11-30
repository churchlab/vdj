#! /usr/bin/env python

import os
import sys
import optparse

import lsf
import seqtools
import vdj
import vdj.pipeline

join = os.path.join

option_parser = optparse.OptionParser()
# option_parser.add_option('-x','--xxx',dest='xxxx',type='int')
(options,args) = option_parser.parse_args()

if len(args) != 1:
    raise ValueError, "Require a single input jobfile"

jobfile = args[0]

# PARAMETER DEFINITION

# process jobfile for input parameters
# defines the following variables:
    # basename               # unique base identifier for data
    # vdj_package_dir        # directory to find scripts.  will be added to shell PATH
    # input_fasta            # the initial fasta data
    # barcode_fasta          # barcode identifiers
    # isotype_fasta          # isotype identifiers
    # analysis_dir           # full path; base directory for final products
    # work_dir               # full path; base directory for intermediate parts and logs
    # min_size               # min size selection
    # max_size               # max size selection
    # packet_size            # packet size for alignment jobs
    # loci                   # the loci to use for VDJ aln
    # raw_vdjxml             # derived file stored in analysis_dir
    # aligned_file           # derived file stored in analysis_dir
    # vj_filtered_file       # derived file stored in analysis_dir
    # size_selected_file     # derived file stored in analysis_dir
    # clustered_file         # derived file stored in analysis_dir
    # parts_dir              # derived intermediate work directories (rel. to work_dir)
    # log_dir                # derived intermediate work directories (rel. to work_dir)
    # partition_dir          # derived intermediate work directories (rel. to work_dir)


execfile(jobfile)

# add script directory to path
os.environ['PATH'] = join(vdj_package_dir,'bin') + ':' + os.environ['PATH']



# make working directories
if not os.path.exists(analysis_dir):                 os.mkdir(analysis_dir,0755)
if not os.path.exists(work_dir):                     os.mkdir(work_dir,0755)
if not os.path.exists(join(work_dir,parts_dir)):     os.mkdir(join(work_dir,parts_dir),0755)
if not os.path.exists(join(work_dir,log_dir)):       os.mkdir(join(work_dir,log_dir),0755)
if not os.path.exists(join(work_dir,partition_dir)): os.mkdir(join(work_dir,partition_dir),0755)

locus_options = ' '.join([' --locus %s' % locus for locus in loci.split()])


# PIPELINE STARTS HERE

# 0. CONVERSION TO VDJXML

sys.stdout.write("Converting FASTA to VDJXML..."); sys.stdout.flush()
inhandle = open(input_fasta,'r')
outhandle = open(join(analysis_dir,raw_vdjxml),'w')
vdj.pipeline.fasta2vdjxml(inhandle,outhandle)
inhandle.close()
outhandle.close()
sys.stdout.write("finished\n"); sys.stdout.flush()

# 1. SIZE SELECTION
sys.stdout.write("Performing size selection on reads..."); sys.stdout.flush()
inhandle = open(join(analysis_dir,raw_vdjxml),'r')
outhandle = open(join(analysis_dir,size_selected_file),'w')
vdj.pipeline.size_select(inhandle,outhandle,min_size,max_size)
inhandle.close()
outhandle.close()
sys.stdout.write("finished\n"); sys.stdout.flush()

# 2. SPLIT INTO PARTS
sys.stdout.write("Splitting input into small parts..."); sys.stdout.flush()
inhandle = open(join(analysis_dir,size_selected_file),'r')
parts = vdj.pipeline.iterator2parts( vdj.parse_VDJXML(inhandle),
                                     join(work_dir,parts_dir,size_selected_file),
                                     packet_size,
                                     prefix='<root>',
                                     suffix='</root>')
sys.stdout.write("finished\n"); sys.stdout.flush()

# BARCODE ID, CODING STRAND, ISOTYPE ID, VDJ CLASSIFICATION via LSF
sys.stdout.write("Setting up LSF command...\n"); sys.stdout.flush()
cmd = 'barcode_id.py --barcodes %s ' % barcode_fasta       # 3. BARCODE IDENTIFICATION
cmd += ' | coding_strand.py' + locus_options               # 4. CODING STRAND
if 'IGH' in loci:                                       # 5. ISOTYPE ID (heavy chain only)
    cmd += ' | isotype_id.py --IGHC %s' % isotype_fasta
cmd += ' | align_vdj.py' + locus_options                   # 6. VDJ CLASSIFICATION

# submit cmd to LSF for each part
sys.stdout.write("Submitting jobs to LSF..."); sys.stdout.flush()
jobIDs = []
logfiles = []
outnames = []
for part in parts:
    partID = part.split('.')[-1]
    partoutname = join(work_dir,parts_dir,basename)+'.aligned.vdjxml.'+partID
    outnames.append(partoutname)
    curr_cmd = 'cat %s | ' + cmd + ' > %s'
    curr_cmd = curr_cmd % (part,partoutname)
    logfile = join(work_dir,log_dir,'alignment.log.')+partID
    jobID = lsf.submit_to_LSF(short_queue,logfile,curr_cmd)
    logfiles.append(logfile)
    jobIDs.append(jobID)
sys.stdout.write("finished\n"); sys.stdout.flush()

sys.stdout.write("Waiting for LSF jobs..."); sys.stdout.flush()
lsf.wait_for_LSF_jobs(jobIDs,logfiles)
sys.stdout.write("finished\n"); sys.stdout.flush()

# 7. CONCAT PARTS
sys.stdout.write("Concatenating pieces..."); sys.stdout.flush()
outhandle = open(join(analysis_dir,aligned_file),'w')
vdj.pipeline.cat_vdjxml(outnames,outhandle)
outhandle.close()
sys.stdout.write("finished\n"); sys.stdout.flush()

# 8. FILTER VJ
sys.stdout.write("Filtering for chains with VJ alignments..."); sys.stdout.flush()
inhandle = open(join(analysis_dir,aligned_file),'r')
outhandle = open(join(analysis_dir,vj_filtered_file),'w')
vdj.pipeline.filter_VJ(inhandle,outhandle)
inhandle.close()
outhandle.close()
sys.stdout.write("finished\n"); sys.stdout.flush()

# 9. PARTITION VJ
sys.stdout.write("Partitioning data by VJ alignment..."); sys.stdout.flush()
inhandle = open(join(analysis_dir,vj_filtered_file),'r')
partition_files = vdj.pipeline.partition_VJ(inhandle,join(work_dir,partition_dir,basename))
inhandle.close()
sys.stdout.write("finished\n"); sys.stdout.flush()

# 10. CLUSTER CDR3s
sys.stdout.write("Setting up clustering jobs...\n"); sys.stdout.flush()
cmd = 'cluster_cdr3.py --cutoff %f --linkage %s' % (clustering_cutoff,clustering_linkage)

# submit LSF
sys.stdout.write("Submitting clustering jobs to LSF..."); sys.stdout.flush()
jobIDs = []
logfiles = []
outnames = []
for partition in partition_files:
    partitionoutname = '.'.join(partition.split('.')[:-1])+'.clustered.vdjxml'
    outnames.append(partitionoutname)
    tag = os.path.splitext(os.path.basename(partition))[0].split('.')[-1]
    curr_cmd = 'cat %s | ' + cmd + ' --tag %s' + ' > %s'
    curr_cmd = curr_cmd % (partition,tag,partitionoutname)
    logfile = join(work_dir,log_dir,'clustering.')+tag+'.log'
    jobID = lsf.submit_to_LSF(long_queue,logfile,curr_cmd)
    logfiles.append(logfile)
    jobIDs.append(jobID)
sys.stdout.write("finished\n"); sys.stdout.flush()

sys.stdout.write("Waiting for LSF jobs..."); sys.stdout.flush()
lsf.wait_for_LSF_jobs(jobIDs,logfiles)
sys.stdout.write("finished\n"); sys.stdout.flush()

# 11. CONCAT PARTITIONS
sys.stdout.write("Concatenating partitions..."); sys.stdout.flush()
outhandle = open(join(analysis_dir,clustered_file),'w')
vdj.pipeline.cat_vdjxml(outnames,outhandle)
outhandle.close()
sys.stdout.write("finished\n"); sys.stdout.flush()
