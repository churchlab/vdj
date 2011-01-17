"""params.py

Define directory and file names that must be manually modified
to point to certain resources.
"""

# HACK.  Figure out better way to refer to this directory
vdj_dir = '/Users/laserson/code/vdj'

# packaged data dir
data_dir = 'data'
IGHV_fasta = 'IGHV.fasta'
IGHD_fasta = 'IGHD.fasta'
IGHJ_fasta = 'IGHJ.fasta'
IGKV_fasta = 'IGKV.fasta'
IGKJ_fasta = 'IGKJ.fasta'
IGLV_fasta = 'IGLV.fasta'
IGLJ_fasta = 'IGLJ.fasta'
TRBV_fasta = 'TRBV.fasta'
TRBD_fasta = 'TRBD.fasta'
TRBJ_fasta = 'TRBJ.fasta'
TRAV_fasta = 'TRAV.fasta'
TRAJ_fasta = 'TRAJ.fasta'
TRDV_fasta = 'TRDV.fasta'
TRDD_fasta = 'TRDD.fasta'
TRDJ_fasta = 'TRDJ.fasta'
TRGV_fasta = 'TRGV.fasta'
TRGJ_fasta = 'TRGJ.fasta'


# The following directory and files will not be packaged with vdj
# but will be computed the first time refseq is imported.  After
# that, it will not be recomputed unless it is forced

# persistent data directory
pickle_dir = 'pickle'

# Relevant LIGM records in pickle format
# If the file exists, refseq will not try to recompute it
# unless it's forced
IGHV_pickle = 'IGHV.pickle'
IGHD_pickle = 'IGHD.pickle'
IGHJ_pickle = 'IGHJ.pickle'
IGKV_pickle = 'IGKV.pickle'
IGKJ_pickle = 'IGKJ.pickle'
IGLV_pickle = 'IGLV.pickle'
IGLJ_pickle = 'IGLJ.pickle'
TRBV_pickle = 'TRBV.pickle'
TRBD_pickle = 'TRBD.pickle'
TRBJ_pickle = 'TRBJ.pickle'
TRAV_pickle = 'TRAV.pickle'
TRAJ_pickle = 'TRAJ.pickle'
TRDV_pickle = 'TRDV.pickle'
TRDD_pickle = 'TRDD.pickle'
TRDJ_pickle = 'TRDJ.pickle'
TRGV_pickle = 'TRGV.pickle'
TRGJ_pickle = 'TRGJ.pickle'


# # full IMGT flatfile database dir
# imgt_dir = '/Users/laserson/research/church/vdj-ome/ref-data/IMGT'
# ligm_filename = 'imgt.dat'

# refdatadir = '/Users/laserson/research/church/vdj-ome/ref-data/IMGT'
# imgtdat = 'imgt.dat'
# imgtfasta = 'imgt.fasta'
# imgtvseq = 'vdj_ref.fasta'
# imgtrefseqfasta = 'imgtrefseq.fasta'
# imgtspecfasta  = 'imgtspec.fasta'
# imgtspecvdjxml = 'imgtspec.vdjxml'
# imgtrefseqFR3endcoords = 'imgtrefseqFR3endcoords.dat'
# imgtrefseqJTRPstartcoords = 'imgtrefseqJTRPstartcoords.dat'
