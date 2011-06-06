import os
import sys
import argparse

argparser = argparse.ArgumentParser(description=None)
argparser.add_argument('inputfile')
argparser.add_argument('--work_dir')
argparser.add_argument('--packet_size',type=int)
argparser.add_argument('--locus',action='append')
argparser.add_argument('--option',dest='xxx',action='store_const',default=5)
args = argparser.parse_args()


# 0. SETUP
join = os.path.join
def log(s): sys.stdout.write(s); sys.stdout.flush()

work_dir = os.path.abspath(args.work_dir)
assert not os.path.exists(args.work_dir)
os.mkdir(work_dir,0755)
os.mkdir(join(work_dir,'parts'),0755)
os.mkdir(join(work_dir,'logs'),0755)



# 1. SPLIT INTO PARTS
log("Splitting input into small parts...")
parts = vdj.pipeline.iterator2parts( vdj.parse_imgt(args.inputfile),
                                     join(work_dir,'parts/input.imgt'),
                                     args.packet_size)
log("finished\n")



# 3-7. BARCODE ID, CODING STRAND, ISOTYPE ID, VDJ CLASSIFICATION, TRANSLATION via LSF
log("Setting up LSF command...\n")
locus_options = ' '.join([' --locus %s' % locus for locus in params['locus']])
cmd = 'align_vdj.py' + locus_options                            # 6. VDJ CLASSIFICATION


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
with open(aligned_file,'w') as outhandle:
    vdj.pipeline.cat_vdjxml(outnames,outhandle)
log("finished\n")
