"""params.py

Define directory and file names that must be manually modified
to point to certain resources.
"""

# directory that contains imgt.dat, imgt.fasta, imgtrefseq.fasta from IMGT
# NOTE: after unzipping IMGT.zip, it bungles imgt.dat and imgt.fasta.
#       These need to be deleted and the .Z version SEPARATELY uncompressed.
#
# NOTE: this should be fixed to not require manual pointing to IMGT

# packaged data dir
data_dir = 'data'
IGHV_filename = 'IGHV.fasta'
IGHD_filename = 'IGHD.fasta'
IGHJ_filename = 'IGHJ.fasta'
IGKV_filename = 'IGKV.fasta'
IGKJ_filename = 'IGKJ.fasta'
IGLV_filename = 'IGLV.fasta'
IGLJ_filename = 'IGLJ.fasta'
TRBV_filename = 'TRBV.fasta'
TRBD_filename = 'TRBD.fasta'
TRBJ_filename = 'TRBJ.fasta'
TRAV_filename = 'TRAV.fasta'
TRAJ_filename = 'TRAJ.fasta'
TRDV_filename = 'TRDV.fasta'
TRDD_filename = 'TRDD.fasta'
TRDJ_filename = 'TRDJ.fasta'
TRGV_filename = 'TRGV.fasta'
TRGJ_filename = 'TRGJ.fasta'

# full IMGT flatfile database dir
imgt_dir = '/Users/laserson/research/church/vdj-ome/ref-data/IMGT'
ligm_filename = 'imgt.dat'

# refdatadir = '/Users/laserson/research/church/vdj-ome/ref-data/IMGT'
# imgtdat = 'imgt.dat'
# imgtfasta = 'imgt.fasta'
# imgtvseq = 'vdj_ref.fasta'
# imgtrefseqfasta = 'imgtrefseq.fasta'
# imgtspecfasta  = 'imgtspec.fasta'
# imgtspecvdjxml = 'imgtspec.vdjxml'
# imgtrefseqFR3endcoords = 'imgtrefseqFR3endcoords.dat'
# imgtrefseqJTRPstartcoords = 'imgtrefseqJTRPstartcoords.dat'
