# vdj_pipeline_sizeselect.py

# gets a baseoutputname, and dumps as a new fasta in baseoutputname.sizeMIN-MAX.fasta
# MIN is 0 by default, MAX is inf

import glob
import optparse
import numpy as np
import seqtools
from Bio import SeqIO

# process command line
cmdlineparse = optparse.OptionParser()
cmdlineparse.add_option('--minreadlen', action='store', type='float',   default=0,      dest='minreadlen')
cmdlineparse.add_option('--maxreadlen', action='store', type='float', default=np.inf, dest='maxreadlen')
cmdlineparse.add_option('--outputname', action='store', type='string', dest='outputname')
options, files = cmdlineparse.parse_args()

minreadlen = options.minreadlen
maxreadlen = options.maxreadlen

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
print "Loaded %i sequences for size selection" % len(seqs)

# filter by size
seqs_size_filtered = [seq for seq in seqs if len(seq) >= minreadlen and len(seq) <= maxreadlen]

# dump in file
currfilename = baseoutputname+'.size'+str(minreadlen)+'-'+str(maxreadlen)+'.fasta'
op = open(currfilename,'w')
SeqIO.write( seqs_size_filtered, op, 'fasta' )
op.close()
print "Processed %i sequences of length %i-%i" % (len(seqs_size_filtered),minreadlen,maxreadlen)