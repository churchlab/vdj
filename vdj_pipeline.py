# vdjpipeline.py

# note: globbed filenames must be put in quotations.  they must be fasta files

import glob
import optparse
import seqtools
from Bio import SeqIO
import Bio.Seq
import os
import subprocess
import numpy as np
import time
import vdj
import alignment
from datetime import datetime

# ========================
# = Process command line =
# ========================

# command line to look like this:
# vdj_pipeline operation relevant_options files_to_operate_on
#
# operations can be one of:
#	full				perform all operations below; takes list of fileglobs
#	initial_import		take fasta file and set up initial ImmuneChain/Repertoire XML file
#						takes list of fileglobs
#	size_select			filter out sequences of a certain size; takes single XML file
#	barcode_id			tag with barcode; takes single XML file
#	positive_strand		convert to positive strand; takes single XML file
#	align				perform VDJ alignment, CDR3 extraction, isotype ID; takes single XML file

cmdlineparse = optparse.OptionParser()

cmdlineparse.add_option('--tag', action='append', type='string',dest='tags')
cmdlineparse.add_option('--metatag', action='append', type='string',dest='metatags')
	# add tags to the Repertoire (metatag) or ImmuneChain objs (tag)



cmdlineparse.add_option('-o','--outputname', action='store', type='string', dest='outputname')
	# filename for output files.  takes it from input file if not given
cmdlineparse.add_option('-s','--sizeselect', action='store', type='float', dest='readlensizes', nargs=2)
	# flag for size selection, takes 2 args: min and max
cmdlineparse.add_option('-b','--barcodesplit', action='store', type='string', dest='barcodefile')
	# note that barcodes must be listed in a file followed by an identifier
	# 454 MIDs are: /Users/laserson/research/church/vdj-ome/seq-data/barcodes/454MID.barcodes
cmdlineparse.add_option('-i','--isotype', action='store', type='string',dest='IGHCfile')
	# contains primer seq (5'->3') and then the isotype identifier on each line

# LSF dispatching.  only applies to positive strand identification or VDJ alignment/CDR3 extraction
cmdlineparse.add_option('--LSFdispatch', action='store', type='string', dest='LSFargs', nargs=2)
	# note that first argument is the queue to submit to and second argument is filename to dump LSF output to
cmdlineparse.add_option('--packetsize', action='store', type='int', default=0, dest='packetsize')



options, args = cmdlineparse.parse_args()

# ====================================
# = Initialization before processing =
# ====================================

# what command am i running?
operation = args[0]

# what files am i operating on?
inputfilelist = []
for filename in args[1:]:
	inputfilelist.extend( glob.glob(filename) )

# what is the base name for the output files?
if options.outputname is not None:
	outputname = options.outputname
else:
	outputname = inputfilelist[0]

# ==================================
# = SWITCH based on chosen command =
# ==================================

if operation == 'initial_import':
	seqs = []
	for filename in inputfilelist:
		seqs.extend(seqtools.getFasta(filename))
	
	chains = [ImmuneChain(descr=seq.description.split()[0],seq=seq.seq,tags=options.tags) for seq in seqs]	
	rep = Repertoire(chains)
	rep.add_metatags(options.metatags)
	rep.add_metatags("Initial_Import : " + timestamp())
	
	op = open(outputname,'w')
	vdj.writeVDJ(rep,op)
	op.close()
	
elif operation == 'size_select':
	if len(inputfilelist) > 1:
		raise Exception, "Too many input files for size_select operation."
	if len(options.readlensizes) != 2:
		raise Exception, "Incorrect number of args for size_select operation."
	
	minreadlen = options.readlensizes[0]
	maxreadlen = options.readlensizes[1]
	
	rep = vdj.readVDJ(inputfilelist[0],mode='Repertoire')
	idxs = [i for (i,chain) in enumerate(rep) if len(chain) >= minreadlen and len(chain) <= maxreadlen]
	
	rep_sizeselected = rep[idxs]
	rep_sizeselected.add_metatags("Size Selection: " +"min " + str(minreadlen)+ "max " + str(maxreadlen) + " : " + timestamp())
	
	op = open(outputname,'w')
	vdj.writeVDJ(rep_sizeselected,op)
	op.close()
	
elif operation == 'barcode_id':
	if len(inputfilelist) > 1:
		raise Exception, "Too many input files for size_select operation."
	
	# load the barcodes from file and process them
	ip = open(options.barcodefile,'r')
	barcodes = {}
	for line in ip:
		bc = line.split()
		barcodes[bc[0]] = bc[1]
	ip.close()
	barcodelen = len(barcodes.keys()[0])
	
	# load data
	rep = vdj.readVDJ(inputfilelist[0],mode='Repertoire')
	rep.add_metatags("Barcode_ID : " + timestamp())
	
	for chain in rep:
		curr_bcID = barcodes.get( chain.seq[:barcodelen], '' )
		chain.seq = chain.seq[barcodelen:]	# remove barcode
		chain.add_tags(curr_bcID)
	
	rep.reprocessTags()
	
	op = open(outputname,'w')
	vdj.writeVDJ(rep,op)
	op.close()
	
elif operation == 'positive_strand':
	if len(inputfilelist) > 1:
		raise Exception, "Too many input files for size_select operation."
	
	# load data
	rep = vdj.readVDJ(inputfilelist[0],mode='Repertoire')
	
	# split into LSF and rerun on pieces through LSF system
	if options.LSFargs is not None:
		# get LSF args
		queue = options.LSFargs[0]
		opfile = options.LSFargs[1]
		packetsize = options.packetsize
		
		processes = []	# list of PIDs that are dispatched through LSF
		parts = []
		nfiles = int(math.ceil(float(len(rep)) / packetsize))
		for i in xrange(nfiles):
			# split into multiple parts; append .partNNN to filename
			currfilename = outputname + '.part%03i' % i
			parts.append(currfilename)
			op = open(currfilename,'w')
			vdj.writeVDJ(rep[i*packetsize:(i+1)*packetsize],op)
			op.close()
			# submit to LSF and save Job ID
			proc = subprocess.Popen( ['bsub','-q'+queue,'-o'+opfile,'python','vdj_pipeline.py','positive_strand',currfilename], stdout=subprocess.PIPE )
			processes.append(proc.stdout.read().split('<')[1].split('>')[0])
		
		# wait for jobs to finish
		waitforLSFjobs(processes,30)
		
		# load parts
		rep = Repertoire()
		for part in parts:
			rep_part = vdj.readVDJ(part,mode='Repertoire')
			rep += rep_part

		op = open(outputname,'w')
		vdj.writeVDJ(rep,op)
		op.close()
		
	# just run the algo on the single file
	else:
		rep.add_metatags("Positive_Strand : " + timestamp())
		
		# compute correct strand
		strands = alignment.chainlist2strandlist(rep)
		if len(strands) != len(rep):
			raise Exception, "Error in positive strand id."
		
		# invert strand
		for (strand,chain) in zip(strands,rep):
			chain.add_tags('positive')
			if strand == -1:
				chain.add_tags('revcomp')
				chain.seq = Bio.Seq.reverse_complement(chain.seq)
		
		rep.reprocessTags()
		
		op = open(outputname,'w')
		vdj.writeVDJ(rep,op)
		op.close()
	
elif operation == 'align':
	if len(inputfilelist) > 1:
		raise Exception, "Too many input files for size_select operation."
	
	# load data
	rep = vdj.readVDJ(inputfilelist[0],mode='Repertoire')
	
	# split into LSF and rerun on pieces through LSF system
	if options.LSFargs is not None:
		# get LSF args
		queue = options.LSFargs[0]
		opfile = options.LSFargs[1]
		packetsize = options.packetsize
		
		processes = []	# list of PIDs that are dispatched through LSF
		parts = []
		nfiles = int(math.ceil(float(len(rep)) / packetsize))
		for i in xrange(nfiles):
			# split into multiple parts; append .partNNN to filename
			currfilename = outputname + '.part%03i' % i
			parts.append(currfilename)
			op = open(currfilename,'w')
			vdj.writeVDJ(rep[i*packetsize:(i+1)*packetsize],op)
			op.close()
			# submit to LSF and save Job ID
			proc = subprocess.Popen( ['bsub','-q'+queue,'-o'+opfile,'python','vdj_pipeline.py','align',currfilename], stdout=subprocess.PIPE )
			processes.append(proc.stdout.read().split('<')[1].split('>')[0])
		
		# wait for jobs to finish
		waitforLSFjobs(processes,30)
		
		# load parts
		rep = Repertoire()
		for part in parts:
			rep_part = vdj.readVDJ(part,mode='Repertoire')
			rep += rep_part

		op = open(outputname,'w')
		vdj.writeVDJ(rep,op)
		op.close()
		
	# just run the algo on the single file
	else:
		rep.add_metatags("VDJ_Alignment : " + timestamp())
		
		# perform alignment
		vdj.alignment.abacus.alignchainlist(rep)
		vdj.clones.extractCDR3(rep)
		
		rep.reprocessTags()
		
		op = open(outputname,'w')
		vdj.writeVDJ(rep,op)
		op.close()
elif operation == 'full':
	pass
	
	

def timestamp():
	return datetime.now.isoformat().split('.')[0]

def waitforLSFjobs(PIDs,interval=30):
	finished == False
	while not finished:
		time.sleep(interval)
		p = subprocess.Popen('bjobs',stdout=subprocess.PIPE)
		status = p.stdout.read().split('\n')
		runningprocesses = [line.split()[0] for line in status if line.split() != [] and line.split()[0] != 'JOBID']
		finished = True
		for pid in PIDs:
			if pid in runningprocesses:
				finished = False
	return
