__init__.py         COMPLETE
alignment.py        COMPLETE
analysis.py         COMPLETE
clustering.py       COMPLETE
legacy.py           DELETED
params.py           COMPLETE
pipeline.py         COMPLETE
plotting.py         DELETED
refseq.py           COMPLETE
refsequtils.py      DELETED
setup.py            COMPLETE

Clean up imports

The process for refseq init is this:

1.  `params.py` lists the fasta files that contain the list of reference
    elements, including their informative headers.
2. `refseq.py` checks whether each type of reference set has already been
    computed. If not, it induces the computation by calling
    `process_IMGT_references` in `refsequtils.py`.
3. `process_IMGT_references` then initializes an object from the fasta file
    and in the case of a V region will parse the full LIGM database to get
    annotations for the CDR3.
4.  

Things to update in `alignment.py`:

*   refV_seqs
*   refJ_seqs
*   refD_seqs
*   refV_offset
*   refJ_offset
*   seqdict2kmers   COMPLETE
*   seqdict2revcompseqdict  COMPLETE
*   seq2kmers   COMPLETE
*   bestalignNW COMPLETE
*   bestalignSW COMPLETE
*   construct_alignment COMPLETE
*   Valign_chain    COMPLETE
*   Jalign_chain    COMPLETE

*   `vdj_aligner_combined`

a = '------agtcacggatcg'
b = '------agtc--ggatcg'


* * *
`__init__.py`, `params.py`, `refseq.py` appear to be in finished shape.
However, I ran into an issue I didn't anticipate: it appears that biopython
loads the annotations dictionary in a very particular way for INSDC formats,
and obliterates any user-defined annotations when writing to an INSDC format.
This is a nonstarter for my file format, so I will have to figure out a
workaround. FIXED.
* * *

import vdj
import vdj.alignment
from Bio import SeqIO
from Bio.Alphabet import generic_dna
iter = SeqIO.parse('heavy_chains.HIV.20101014.aligned.reannotated.fasta','fasta',generic_dna)
record = iter.next()
chain = vdj.ImmuneChain(record)
aligner = vdj.alignment.igh_aligner()
aligner.Valign_chain(chain)
reload(vdj)
reload(vdj.alignment)
aligner = vdj.alignment.igh_aligner()
aligner.Valign_chain(chain)

