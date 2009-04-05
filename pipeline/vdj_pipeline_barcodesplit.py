# vdj_pipeline_barcodesplit.py

import glob
import optparse
import seqtools
from Bio import SeqIO

# process command line options
cmdlineparse = optparse.OptionParser()
cmdlineparse.add_option('--barcodefile', action='store', type='string', dest='barcodefile')
	# note that barcodes must be listed in a file followed by an identifier
	# 454 MIDs are: /Users/laserson/research/church/vdj-ome/seq-data/barcodes/454MID.barcodes
cmdlineparse.add_option('--outputname', action='store', type='string', dest='outputname')
options, files = cmdlineparse.parse_args()

# what files am i operating on?
inputfilelist = []
for filename in files:
	inputfilelist.extend( glob.glob(filename) )

# what is the base name for the output files?
if options.outputname is not None:
	baseoutputname = options.outputname
else:
	baseoutputname = '.'.join( inputfilelist[0].split('.')[:-1] )

# load the barcodes from file and process them
ip = open(options.barcodefile,'r')
barcodes = {}
for line in ip:
	bc = line.split()
	barcodes[bc[0]] = bc[1]
ip.close()
barcodelen = len(barcodes.keys()[0])

# load data
seqs = []
for filename in inputfilelist:
	seqs.extend(seqtools.getFasta(filename))
print "Loaded %i sequences for barcode split" % len(seqs)

# split into barcodes
# init:
seqs_split_barcodes = {}	# dict: key is barcode id, and value is list of seqs
seqs_split_barcodes[''] = []	# for reads that fail to match a barcode
for (bcseq,bcID) in barcodes.iteritems():
	seqs_split_barcodes[bcID] = []
# process sequences:
for seq in seqs:
	curr_bcID = barcodes.get( seqtools.seqString(seq)[:barcodelen], '' )
	seq.seq = seq.seq[barcodelen:]	# remove barcode
	seq.description = seq.description.split()[0]
	seqs_split_barcodes[curr_bcID].append(seq)

# dump into files
print "Split sequences into barcodes:"
for (bcID,seqs) in seqs_split_barcodes.iteritems():
	if len(seqs) > 0:
		op = open(baseoutputname+'.barcode_'+bcID+'.fasta','w')
		SeqIO.write(seqs,op,'fasta')
		op.close()
		print bcID + '\t' + str(len(seqs))