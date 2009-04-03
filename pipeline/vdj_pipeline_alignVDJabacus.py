# vdj_pipeline_alignVDJabacus.py

# note: best to put explicit outputname

import glob
import optparse
import seqtools
from Bio import SeqIO
import vdj

# process command line
cmdlineparse = optparse.OptionParser()
cmdlineparse.add_option('--outputname', action='store', type='string', dest='outputname')
cmdlineparse.add_option('--tags1', action='store', type='string',dest='tags',nargs=1)
cmdlineparse.add_option('--tags2', action='store', type='string',dest='tags',nargs=2)
cmdlineparse.add_option('--tags3', action='store', type='string',dest='tags',nargs=3)
cmdlineparse.add_option('--tags4', action='store', type='string',dest='tags',nargs=4)
cmdlineparse.add_option('--tags5', action='store', type='string',dest='tags',nargs=5)
options, files = cmdlineparse.parse_args()

# collect tags
tags = set([])
tags.update(options.tags1)
tags.update(options.tags2)
tags.update(options.tags3)
tags.update(options.tags4)
tags.update(options.tags5)

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
print "Loaded %i sequences for VDJ alignment" % len(seqs)

# perform alignment
chains = vdj.alignment.abacus.alignseqlist(seqs)

rep = vdj.Repertoire(chains)
rep.set_tags(tags)	# add tags

# dump in file
currfilename = baseoutputname+'.vdj'+'.chains'
op = open(currfilename,'w')
vdj.writeImmuneChains(rep,op)
op.close()
print "Aligned %i sequences" % len(rep)