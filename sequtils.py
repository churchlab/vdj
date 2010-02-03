import string

import Bio.SeqIO

revcomp_from = 'ACGTRYSWKMBDHVNacgtryswkmbdhvn'
revcomp_to   = 'TGCAYRSWMKVHDBNtgcayrswmkvhdbn'
revcomp_table = string.maketrans(revcomp_from,revcomp_to)
def reverse_complement(seq):
    return seq.translate(revcomp_table)[::-1]

def load_fasta(filename):
    # makes no guarantee as to the case of the sequence
    ip = open(filename,'r')
    seqs = list(Bio.SeqIO.parse(ip,'fasta'))
    ip.close()
    return seqs

# for generating 'safe' filenames from identifiers
cleanup_table = string.maketrans('/*|','___')
def cleanup_id(identifier):
    return identifier.translate(cleanup_table)
