import sys
import string
from cStringIO import StringIO
import cPickle as pickle
import urllib2

import warnings # this is because IMGT has tons of errors

import ClientForm
from Bio import SeqIO

import seqtools
import params

identity = string.maketrans('','')

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
    partial will take values of "partial in 5'", "partial in 3'", or "partial in 5' and 3'"
    depending where the del is.
    """
    
    def __init__(self,**kw):
        def kw_init(attrib):
            if kw.has_key(attrib):
                self.__setattr__(attrib,kw[attrib])
        
        kw_init('accession')
        kw_init('seq')
        kw_init('gapped_seq')
        kw_init('description')
        kw_init('locus')
        kw_init('gene')
        kw_init('allele')
        kw_init('species')
        kw_init('functional')
        kw_init('imgt_label')
        kw_init('accession_coords')
        kw_init('length')
        kw_init('frame')
        kw_init('partial')
    
    def init_from_imgt(self,fasta_header,seq):
        """Initialize object from IMGT fasta header and seq"""
        data = fasta_header.lstrip('>').rstrip().split('|')
        self.accession = data[0]
        self.gapped_seq = seq
        # self.seq = seq.translate(None,'.').upper()  # remove periods (gaps)   # python >2.5
        self.seq = seq.translate(identity,'.').upper()
        self.description = data[1]
        self.locus = self.description[0:4]
        self.gene = self.description.split('*')[0]
        self.allele = self.description
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
        self.partial = data[13]
    
    def pull_LIGM_record(self):
        """Get SeqRecord object for LIGM record from IMGT server"""
        
        # NOTE: this can potentially be significantly simplified by accessing the URL
        # interface to LIGM, through:
        # http://imgt.cines.fr/cgi-bin/IMGTlect.jv?query=5+numacc
        # where numacc is the accession number
        
        request = urllib2.Request('http://imgt.cines.fr/cgi-bin/IMGTlect.jv?livret=0')
        # LIGM page
        response = urllib2.urlopen(request)
        forms = ClientForm.ParseResponse(response,
                                         form_parser_class=ClientForm.XHTMLCompatibleFormParser,
                                         backwards_compat=False)
        form = forms[1]
        form['l01p01c02'] = self.accession
        request2 = form.click()
        # data format page
        response2 = urllib2.urlopen(request2)
        forms2 = ClientForm.ParseResponse(response2,
                                         form_parser_class=ClientForm.XHTMLCompatibleFormParser,
                                         backwards_compat=False)
        form2 = forms2[0]
        assert( form2.controls[8].attrs['value'] == '2 IMGT flat-file' )
        form2.controls[8].id = 'flatfile'
        request3 = form2.click(id='flatfile')
        # LIGM record results
        response3 = urllib2.urlopen(request3)
    
        # ghetto parse of the results.  the text of the LIGM record is in <pre>...</pre> tags
        rawdata1 = response3.read()
        rawdata2 = rawdata1.split('<pre>')[1].split('</pre>')[0].lstrip()
        rawdata3 = StringIO(rawdata2)
        self.record = SeqIO.read(rawdata3,'imgt')

class VReferenceEntry(ReferenceEntry):
    def __init__(self,**kw):
        ReferenceEntry.__init__(self,**kw)
    
    def set_CDR3_boundary(self):    # FR3 end
        """Get coord of end of FR3 from IMGT LIGM database."""
        # some records can have multiple references in them
        target_allele = self.allele
        feature_iter = self.record.features.__iter__()
        v_gene = seqtools.advance_to_features(feature_iter,['V-REGION','V-GENE'])
        while v_gene.qualifiers.get('allele',[''])[0] != target_allele:  # advance to the target gene
            v_gene = seqtools.advance_to_features(feature_iter,['V-REGION','V-GENE'])
        conserved_cys = seqtools.advance_to_feature(feature_iter,'2nd-CYS')
        
        # note: biopython features already use pythonic indexing
        self.cdr3_boundary = conserved_cys.location.start.position

class JReferenceEntry(ReferenceEntry):
    def __init__(self,**kw):
        ReferenceEntry.__init__(self,**kw)
        
    def set_CDR3_boundary(self):    # FR4 start
        """Get coord of start of FR4 from IMGT LIGM database."""
        # some records can have multiple references in them
        target_allele = self.allele
        feature_iter = self.record.features.__iter__()
        j_gene = seqtools.advance_to_features(feature_iter,['J-REGION','J-GENE'])
        while j_gene.qualifiers.get('allele',[''])[0] != target_allele:  # advance to the target gene
            j_gene = seqtools.advance_to_features(feature_iter,['J-REGION','J-GENE'])
        # note, there can be a conserved TRP or PHE
        conserved_trp = seqtools.advance_to_features(feature_iter,['J-TRP','J-PHE'])
        
        # note: biopython features already use pythonic indexing
        self.cdr3_boundary = conserved_trp.location.end.position


# ================
# = Parsing IMGT =
# ================

# import pdb

def process_IMGT_references(ref_entry_cls,fasta_infilename,pickle_outfilename,verbose=False):
    """Load references from the IMGT/V-QUEST fasta file refs
    
    e.g., IGHV.fasta, present in the data directory
    ref_entry_cls is the class object for the reference type,
        e.g., VReferenceEntry, JReferenceEntry, etc.
    """
    
    # pdb.set_trace()
    
    references = []
    ip = open(fasta_infilename,'r')
    for record in SeqIO.parse(ip,'fasta'):
        curr_reference = ref_entry_cls()
        curr_header = record.description
        curr_seq = record.seq.tostring()
        
        # Potential problems with FASTA headers
        try:
            curr_reference.init_from_imgt(curr_header,curr_seq)
        except ValueError:
            warnings.warn("Invalid header: %s" % curr_header)
            continue
        
        # I don't want to deal with partial seqs right now
        if 'partial' in curr_reference.partial:
            continue
        
        curr_reference.pull_LIGM_record()
        # Potential problems with finding annotated CDR3 boundary
        try:
            curr_reference.set_CDR3_boundary()
        except AttributeError:
            # ReferenceEntry has no set_CDR3_boundary. Used for D segments
            pass
        except ValueError, e:
            warnings.warn("Failed to find CDR3 boundary in %s. Skipping..." % curr_reference.allele)
            continue
        
        references.append(curr_reference)
        if verbose: print "Finished processing %s" % curr_reference.allele
        sys.stdout.flush()
    ip.close()
    
    op = open(pickle_outfilename,'w')
    pickle.dump(references,op,protocol=2)
    op.close()
    
    return references
