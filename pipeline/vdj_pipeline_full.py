# vdjpipeline.py

# note: globbed filenames must be put in quotations.  they must be fasta files

import glob
import optparse
import seqtools
from Bio import SeqIO
import numpy as np
import os

# size select
# split barcodes
# positive strand

# each of the preceding gives a fasta file with an additional keyword before .fasta
# the next one switches from .fasta to .chains, so we need to give an explicit name
# and supply tags as well to incorporate into the ImmuneChain objects

# align VDJ
# get isotype
# extract CDR3

# finally, recombine all sequences together appropriately and delete intermediate parts
