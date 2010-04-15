class ReferenceEntry(object):
    """Data structure to hold a reference sequence from IMGT/GENE-DB or
    IMGT/V-QUEST.
    
    Some of the attributes are computed, and some are taken from IMGT.
    
    accession_coords uses converted python numbering
    gapped_seq includes IMGT gaps that they provide
    frame uses python coords: 0 means already in frame.
                              1 means skip one nt
                              2 means skip two nts
    partial will take values of "5'" or "3'", depending where the del is.
    """
    
    def __init__(self,
                 seq=None,
                 gapped_seq=None,
                 accession=None,
                 description=None,
                 locus=None,
                 gene=None,
                 allele=None,
                 species=None,
                 functional=None,
                 imgt_label=None,
                 accession_coords=None,
                 length=None,
                 frame=None,
                 partial=None):
        self.accession = accession
        self.seq = seq
        self.gapped_seq = gapped_seq
        self.description = description
        self.locus = locus
        self.gene = gene
        self.allele = allele
        self.species = species
        self.functional = functional
        self.imgt_label = imgt_label
        self.accession_coords = accession_coords
        self.length = length
        self.frame = frame
        self.partial = partial
    
    def init_from_imgt(self,fasta_header,seq):
        """Initialize object from IMGT fasta header and seq"""
        data = fasta_header.lstrip('>').split('|')
        self.accession = data[0]
        self.gapped_seq = seq
        self.seq = seq.translate(None,'.')  # remove periods (gaps)
        self.description = data[1]
        self.locus = self.description[0:4]
        self.gene = self.description.split('*')[0]
        self.allele = self.description.split('*')[1]
        self.species = data[2]
        self.functional = data[3]
        self.imgt_label = data[4]
        raw_coords_start = int(data[5].split('.')[0])
        raw_coords_end   = int(data[5].split('.')[-1])
        coords = (raw_coords_start - 1,raw_coords_end)    # note the conversion to python coord system
        self.accession_coords = coords
        self.length = len(self.seq)
        if self.length != int(data[6].split()[0]):
            raise ValueError, "Lengths are inconsistent: %s" % fasta_header
        self.frame = int(data[7]) - 1   # note change to python numbering (0-based)
        self.partial = data[13].split()[-1]


# TODO:
# Subclass the ReferenceEntry for V and J, which contain addl features of interest
# such as conserved internal positions (e.g., 2nd-CYS)