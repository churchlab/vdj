__init__.py         COMPLETE
alignment.py        
analysis.py         COMPLETE
clustering.py       COMPLETE
legacy.py           DELETED
params.py           COMPLETE
pipeline.py         
plotting.py         DELETED
refseq.py           COMPLETE
refsequtils.py      
setup.py            COMPLETE

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
*   seqdict2kmers
*   seqdict2revcompseqdict
*   seq2kmers
*   bestalignNW
