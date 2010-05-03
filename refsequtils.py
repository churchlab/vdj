# ===================
# = Data structures =
# ===================

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
    
    def __init__(self,**kw):
        # self.seq=None
        # self.gapped_seq=None
        # self.accession=None
        # self.description=None
        # self.locus=None
        # self.gene=None
        # self.allele=None
        # self.species=None
        # self.functional=None
        # self.imgt_label=None
        # self.accession_coords=None
        # self.length=None
        # self.frame=None
        # self.partial=None
        
        if kw.has_key('accession'): self.accession = accession
        if kw.has_key('seq'): self.seq = seq
        if kw.has_key('gapped_seq'): self.gapped_seq = gapped_seq
        if kw.has_key('description'): self.description = description
        if kw.has_key('locus'): self.locus = locus
        if kw.has_key('gene'): self.gene = gene
        if kw.has_key('allele'): self.allele = allele
        if kw.has_key('species'): self.species = species
        if kw.has_key('functional'): self.functional = functional
        if kw.has_key('imgt_label'): self.imgt_label = imgt_label
        if kw.has_key('accession_coords'): self.accession_coords = accession_coords
        if kw.has_key('length'): self.length = length
        if kw.has_key('frame'): self.frame = frame
        if kw.has_key('partial'): self.partial = partial
    
    def init_from_imgt(self,fasta_header,seq):
        """Initialize object from IMGT fasta header and seq"""
        data = fasta_header.lstrip('>').split('|')
        self.accession = data[0]
        self.gapped_seq = seq
        self.seq = seq.translate(None,'.').upper()  # remove periods (gaps)
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


class VReferenceEntry(ReferenceEntry):
    def __init__(self,**kw):
        self.ReferenceEntry.__init__(**kw)


class JReferenceEntry(ReferenceEntry):
    def __init__(self,**kw):
        self.ReferenceEntry.__init__(**kw)



# ================
# = Parsing IMGT =
# ================


import os
import simplejson as json

from Bio import SeqIO

def load_references(locus,data_dir,filename):
    references = []
    ip = open(os.path.join(data_dir,filename),'r')
    for record in SeqIO.parse(ip,'fasta'):
        if locus[-1] == 'V':
            curr_reference = VReferenceEntry()
            curr_header = record.description
            curr_seq = record.seq.tostring()
            curr_reference.init_from_imgt(curr_header,curr_seq)
            references.append(curr_reference)
        elif locus[-1] == 'J':
            curr_reference = JReferenceEntry()
            curr_header = record.description
            curr_seq = record.seq.tostring()
            curr_reference.init_from_imgt(curr_header,curr_seq)
            references.append(curr_reference)
        else:
            curr_reference = ReferenceEntry()
            curr_header = record.description
            curr_seq = record.seq.tostring()
            curr_reference.init_from_imgt(curr_header,curr_seq)
            references.append(curr_reference)
    ip.close()


def get_LIGM_records(accessions,data_dir,imgt_dir,ligm_file,output_file):
    records = {}
    
    ip = open(os.path.join(imgt_dir,ligm_file),'r')
    for record in SeqIO.parse(ip,'embl'):
        if name in accessions:
            records[name] = record
    ip.close()
    
    op = open(output_file,'w')
    json.dump(records,op)
    op.close()
    
    return records


def get_FR3_IMGT_end(v_ref_elt,ligm_record):
    












