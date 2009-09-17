import string

import Bio.SeqIO

revcomp_from = 'ACGTRYSWKMBDHVN'
revcomp_to   = 'TGCAYRSWMKVHDBN'
revcomp_table = string.maketrans(revcomp_from,revcomp_to)
def reverse_complement(seq):
    return seq.translate(revcomp_table)

def load_fasta(filename):
    ip = open(filename,'r')
    seqs = list(Bio.SeqIO.parse(ip,'fasta'))
    ip.close()
    return seqs