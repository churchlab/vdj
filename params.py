"""params.py

Define directory and file names that must be manually modified
to point to certain resources.
"""

# directory that contains imgt.dat, imgt.fasta, imgtrefseq.fasta from IMGT
# NOTE: after unzipping IMGT.zip, it bungles imgt.dat and imgt.fasta.
#       These need to be deleted and the .Z version SEPARATELY uncompressed.
refdatadir = '/Users/laserson/research/church/vdj-ome/ref-data/IMGT'
imgtdat = 'imgt.dat'
imgtfasta = 'imgt.fasta'
imgtvseq = 'vdj_ref.fasta'
imgtrefseqfasta = 'imgtrefseq.fasta'
imgtspecfasta  = 'imgtspec.fasta'
imgtspecvdjxml = 'imgtspec.vdjxml'
imgtrefseqFR3endcoords = 'imgtrefseqFR3endcoords.dat'
imgtrefseqJTRPstartcoords = 'imgtrefseqJTRPstartcoords.dat'
