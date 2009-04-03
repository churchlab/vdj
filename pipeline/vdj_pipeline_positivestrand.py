# vdj_pipeline_positivestrand.py

import glob
import optparse
import seqtools
from Bio import SeqIO
import vdj

# process command line
cmdlineparse = optparse.OptionParser()
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

# load data
seqs = []
for filename in inputfilelist:
	seqs.extend(seqtools.getFasta(filename))
print "Loaded %i sequences for positive strand identification" % len(seqs)

# ensure only seq.description is nonempty
for seq in seqs:
	seq.id = ''
	seq.name = ''

# find positive strand and revcomp (changes description too)
seqs_posstrand = vdj.alignment.seqlist2positiveseqlist(seqs)

# dump in file
currfilename = baseoutputname+'.posstrand'+'.fasta'
op = open(currfilename,'w')
SeqIO.write( seqs_posstrand, op, 'fasta' )
op.close()
print "Processed %i sequences into positive strand" % len(seqs_posstrand)