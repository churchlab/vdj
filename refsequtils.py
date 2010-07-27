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
        # self.seq                = kw.get('seq',None)
        # self.gapped_seq         = kw.get('gapped_seq',None)
        # self.accession          = kw.get('accession',None)
        # self.description        = kw.get('description',None)
        # self.locus              = kw.get('locus',None)
        # self.gene               = kw.get('gene',None)
        # self.allele             = kw.get('allele',None)
        # self.species            = kw.get('species',None)
        # self.functional         = kw.get('functional',None)
        # self.imgt_label         = kw.get('imgt_label',None)
        # self.accession_coords   = kw.get('accession_coords',None)
        # self.length             = kw.get('length',None)
        # self.frame              = kw.get('frame',None)
        # self.partial            = kw.get('partial',None)
        
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
        data = fasta_header.lstrip('>').rstrip().split('|')
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

import seqtools

def load_IMGT_references(ref_entry_cls,filename):
    """Load references from the IMGT/V-QUEST fasta file refs
    
    e.g., IGHV.fasta, present in the data directory
    ref_entry_cls is the class object for the reference type,
        e.g., VReferenceEntry, JReferenceEntry, etc.
    """
    references = []
    ip = open(filename,'r')
    for record in SeqIO.parse(ip,'fasta'):
        curr_reference = ref_entry_cls()
        curr_header = record.description
        curr_seq = record.seq.tostring()
        curr_reference.init_from_imgt(curr_header,curr_seq)
        references.append(curr_reference)
    ip.close()
    return references

def get_LIGM_records(accessions,ligm_file):
    """Pull set of accessions from full LIGM file imgt.dat"""
    records = {}
    ligm = SeqIO.index(ligm_file,'imgt')
    for name in accessions:
        records[name] = ligm[name]
    return records

def store_LIGM_records(accessions,ligm_file,output_file):
    """Write set of accessions from LIGM as json"""
    records = get_LIGM_records(accessions,ligm_file)
    op = open(output_file,'w')
    json.dump(records,op)
    op.close()

def set_FR3_IMGT_end(v_ref_elt,ligm_record):
    """Get coord of end of FR3 from IMGT LIGM database.
    
    v_ref_elt is a VReferenceEntry object, and ligm_record is its
    corresponding LIGM entry in imgt.dat (biopython SeqRecord obj)
    """
    conserved_cys = seqtools.get_features(ligm_record.features,'2nd-CYS')
    
    # error checking
    if len(conserved_cys) == 0:
        raise ValueError, "could not find 2nd-CYS in %s" % v_ref_elt.accession
    elif len(conserved_cys) > 1:
        raise ValueError, "found multiple 2nd-CYS in %s" % v_ref_elt.accession
    conserved_cys = conserved_cys[0]
    
    # note: biopython features already use pythonic indexing
    v_ref_elt.fr3_end = conserved_cys.location.start.position
    
    return

def get_J_TRP_start(j_ref_elt,ligm_record):
    """Get coord of start of FR4 from IMGT LIGM database.
    
    j_ref_elt is a JReferenceEntry object, and ligm_record is its
    corresponding LIGM entry in imgt.dat (biopython SeqRecord obj)
    
    note: complication because one record can have multiple genes
    """
    # find correct feature
    feature_iter = ligm_record.features.__iter__()
    conserved_trp = seqtools.get_features(ligm_record.features,'2nd-CYS')












