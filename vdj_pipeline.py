# vdjpipeline.py

# note: globbed filenames must be put in quotations.  they must be fasta files

import glob
import optparse
import os

import vdj
import vdj.LSF

# ========================
# = Process command line =
# ========================

# command line to look like this:
# vdj_pipeline operation relevant_options files_to_operate_on
#
# operations can be one of:
#	full				perform all operations below except clustering; takes list of fileglobs
#	initial_import		take fasta file and set up initial ImmuneChain/Repertoire XML file
#						takes list of fileglobs
#	size_select			filter out sequences of a certain size; takes single XML file
#	barcode_id			tag with barcode; takes single XML file
#	isotype_id			tag with isotype; takes single XML file
#	positive_strand		convert to positive strand; takes single XML file
#	align_rep			perform VDJ alignment, CDR3 extraction, isotype ID; takes single XML file
#	cluster_rep			cluster repertoire

cmdlineparse = optparse.OptionParser()
cmdlineparse.add_option('--tag', action='append', type='string',dest='tags',default=[])
cmdlineparse.add_option('--metatag', action='append', type='string',dest='metatags',default=[])
	# add tags to the Repertoire (metatag) or ImmuneChain objs (tag)
cmdlineparse.add_option('-o','--outputname', action='store', type='string', dest='outputname')
	# filename for output files.  takes it from input file if not given
cmdlineparse.add_option('-s','--sizeselect', action='store', type='float', dest='readlensizes', nargs=2)
	# flag for size selection, takes 2 args: min and max
cmdlineparse.add_option('-b','--barcodesplit', action='store', type='string', dest='barcodefile')
	# note that barcodes must be listed in a file followed by an identifier
	# 454 MIDs are: /Users/laserson/research/church/vdj-ome/seq-data/barcodes/454MID.barcodes
cmdlineparse.add_option('-i','--isotype', action='store', type='string', dest='IGHCfile')
	# contains primer seq (5'->3') and then the isotype identifier on each line
	# 20080924 primers are: /Users/laserson/research/church/vdj-ome/primer-sets/isotypeIDs/20080924.IGHC.ID
cmdlineparse.add_option('--cutoff', action='store', type='float', dest='cutoff')
	# cutoff value for clipping clusters off linkage tree; good default it 4.5, but must be explicitly given
cmdlineparse.add_option('--clustertag', action=store, type='string', dest='clustertag', default='')
	# contains a tag that will be added to each clustertag in the clustering operation
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

# ====================================
# = SWITCH based on chosen operation =
# ====================================

if operation == 'initial_import':
	vdj.initial_import(inputfilelist,outputname,options.metatags,options.tags,tag_rep=True)
elif operation == 'size_select':
	if len(inputfilelist) > 1:
		raise Exception, "Too many input files for size_select operation."
	rep = vdj.fastreadVDJ(inputfilelist[0],mode='Repertoire')
	rep = vdj.size_select(rep,options.readlensizes,tag_rep=True)
	vdj.writeVDJ(rep,outputname)
elif operation == 'barcode_id':
	if len(inputfilelist) > 1:
		raise Exception, "Too many input files for size_select operation."
	rep = vdj.fastreadVDJ(inputfilelist[0],mode='Repertoire')
	rep = vdj.barcode_id(rep,options.barcodefile,tag_rep=True)
	vdj.writeVDJ(rep,outputname)
elif operation == 'isotype_id':
	if len(inputfilelist) > 1:
		raise Exception, "Too many input files for size_select operation."
	rep = vdj.fastreadVDJ(inputfilelist[0],mode='Repertoire')
	rep = vdj.isotype_id(rep,options.IGHCfile,tag_rep=True)
	vdj.writeVDJ(rep,outputname)
elif operation == 'positive_strand':
	if len(inputfilelist) > 1:
		raise Exception, "Too many input files for size_select operation."
	if options.LSFargs is not None: # if dispatching to LSF
		rep = vdj.fastreadVDJ(inputfilelist[0],mode='Repertoire')
		parts = vdj.LSF.split_into_parts(rep,outputname,options.packetsize)
		scriptname = vdj.LSF.generate_script(operation)
		processes = vdj.LSF.submit_to_LSF(options.LSFargs[0],options.LSFargs[1],scriptname,parts)
		vdj.LSF.waitforLSFjobs(processes,30)
		if os.path.exists(scriptname): os.remove(scriptname)
		rep = vdj.LSF.load_parts(parts)
		rep.add_metatags("Positive_Strand : " + vdj.timestamp())
		vdj.writeVDJ(rep,outputname)
		for part in parts:
			if os.path.exists(part): os.remove(part)
	else: # if just running on one file
		rep = vdj.fastreadVDJ(inputfilelist[0],mode='Repertoire')		
		rep = vdj.positive_strand(rep,tag_rep=True)
		vdj.writeVDJ(rep,outputname)
elif operation == 'align_rep':
	if len(inputfilelist) > 1:
		raise Exception, "Too many input files for size_select operation."
	if options.LSFargs is not None: # if dispatching to LSF
		rep = vdj.fastreadVDJ(inputfilelist[0],mode='Repertoire')
		parts = vdj.LSF.split_into_parts(rep,outputname,options.packetsize)
		scriptname = vdj.LSF.generate_script(operation)
		processes = vdj.LSF.submit_to_LSF(options.LSFargs[0],options.LSFargs[1],scriptname,parts)
		vdj.LSF.waitforLSFjobs(processes,30)
		if os.path.exists(scriptname): os.remove(scriptname)
		rep = vdj.LSF.load_parts(parts)
		rep.add_metatags("VDJ_Alignment : " + vdj.timestamp())
		vdj.writeVDJ(rep,outputname)
		for part in parts:
			if os.path.exists(part): os.remove(part)
	else:
		rep = vdj.fastreadVDJ(inputfilelist[0],mode='Repertoire')
		rep = vdj.align_rep(rep,tag_rep=True)
		vdj.writeVDJ(rep,outputname)
elif operation == 'cluster_rep':
	if len(inputfilelist) > 1:
		raise Exception, "Too many input files for size_select operation."
	if options.LSFargs is not None: # if dispatching to LSF
		rep = vdj.fastreadVDJ(inputfilelist[0],mode='Repertoire')
		(inparts,vjcombos) = vdj.split_into_good_VJCDR3s(rep,outputname,verbose=True)
		scriptname = vdj.LSF.generate_script(operation,[options.cutoff,options.clustertag])
		processes = vdj.LSF.submit_to_LSF(options.LSFargs[0],options.LSFargs[1],scriptname,parts)
		vdj.LSF.waitforLSFjobs(processes,30)
		if os.path.exists(scriptname): os.remove(scriptname)
		rep = vdj.LSF.load_parts(parts)
		rep.add_metatags("Clustering|" + options.clustertag + "levenshtein|single_linkage|cutoff="+str(options.cutoff)+"|"+timestamp())
		vdj.writeVDJ(rep,outputname)
		for part in parts:
			if os.path.exists(part): os.remove(part)
	else: # if just running on one file serially
		rep = vdj.fastreadVDJ(inputfilelist[0],mode='Repertoire')
		rep = vdj.clusterRepertoire(rep,options.cutoff,tag_rep=True,tag_chains=True,tag=options.clustertag)
		vdj.writeVDJ(rep,outputname)
elif operation == 'full':
	rep = vdj.initial_import(inputfilelist,outputname,options.metatags,options.tags,tag_rep=True)
	#DEBUG
	#print "post initial import: " + str(type(rep))
	#print "post initial import: len of rep: " + str(len(rep))
	#print rep
	rep = vdj.size_select(rep,options.readlensizes,tag_rep=True)
	rep = vdj.barcode_id(rep,options.barcodefile,tag_rep=True)
	if options.LSFargs is not None: # if dispatching to LSF
		parts = vdj.LSF.split_into_parts(rep,outputname,options.packetsize)
		del rep
		scriptname1 = vdj.LSF.generate_script('positive_strand')
		processes1 = vdj.LSF.submit_to_LSF(options.LSFargs[0],options.LSFargs[1],scriptname1,parts)
		vdj.LSF.waitforLSFjobs(processes1,30)
		if os.path.exists(scriptname1): os.remove(scriptname1)
		scriptname2 = vdj.LSF.generate_script('align_rep')
		processes2 = vdj.LSF.submit_to_LSF(options.LSFargs[0],options.LSFargs[1],scriptname2,parts)
		vdj.LSF.waitforLSFjobs(processes2,30)
		if os.path.exists(scriptname2): os.remove(scriptname2)
		rep = vdj.LSF.load_parts(parts)
		rep.add_metatags("Positive_Strand : " + vdj.timestamp())
		rep.add_metatags("VDJ_Alignment : " + vdj.timestamp())
		rep = vdj.isotype_id(rep,options.IGHCfile,tag_rep=True)
		vdj.writeVDJ(rep,outputname)
		for part in parts:
			if os.path.exists(part): os.remove(part)
	else: # if running on a single file
		rep = vdj.positive_strand(rep,tag_rep=True)
		rep = vdj.isotype_id(rep,options.IGHCfile,tag_rep=True)
		rep = vdj.align_rep(rep,tag_rep=True)
		vdj.writeVDJ(rep,outputname)