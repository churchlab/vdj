Lots of scripts for performing vdj operations.


size_select.py              Size select reads
barcode_id.py               Annotate barcode onto sequences
coding_strand.py            Convert chains to coding sequence
isotype_id.py               Annote isotypes
align_vdj.py                Perform vdj classification
filter_VJ.py                Select only chains with V and J alignments
cluster_cdr3.py             Perform hierarchical clustering of ImmuneChains using their junctions
partition_VJ.py             Partitions vdjxml into files by VJ combo


Older generation:

update_vdjxml.py            Update vdjxml from older version to newer version
fasta2vdjxml.py             Convert fasta file to vdjxml. Takes first white-space delim field for descr
vdjxml2parts.py             Split vdjxml file into parts
cat_vdjxml.py               cat operation on vdjxml files (handles root elements)
cluster_split_VJ.py
cluster_split_VJ_LSF.py
filter_tags_and.py
filter_tags_not.py
filter_tags_or.py
split_on_tags.py
tag_chains.py
vdj_full_pipeline_LSF.py
vdjxml2clone_counts.py
vdjxml2fasta.py
