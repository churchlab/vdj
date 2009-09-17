import string

revcomp_from = 'ACGTRYSWKMBDHVN'
revcomp_to   = 'TGCAYRSWMKVHDBN'
revcomp_table = string.maketrans(revcomp_from,revcomp_to)
def reverse_complement(seq):
    return seq.translate(revcomp_table)