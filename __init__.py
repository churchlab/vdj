import types

from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import NucleotideAlphabet
from Bio import SeqIO

from SeqRecordLite import SeqRecordLite

# select which IMGT reference set to load; human by default
organism = 'human'
# organism = 'mouse'

# ===================
# = DATA STRUCTURES =
# ===================

class ImmuneChain(SeqRecordLite):
    """Data structure to represent an immune chain.
    
    The underlying data structure is actually a biopython `SeqRecord` object.
    This class wraps it in a way that maintains the simple-to-use interface to
    get at some common annotations. 
    
    Sequences must be nucleotides, not protein.
    """
    
    def __init__(self,**kw):
        """Initialize ImmuneChain
        
        seq is 5'->3'
        """
        # first we define our underlying SeqRecord object
        if 'record' in kw:  # init with SeqRecord object
            SeqRecordLite.__init__(self,kw['record'])
        else:   # otherwise, initialize manually
            if 'seq' in kw:
                seq = Seq(data=kw['seq'],alphabet=NucleotideAlphabet())
                record = SeqRecord(seq=seq,id='',name='',description='')
                SeqRecordLite.__init__(self,record)
            else:
                SeqRecordLite.__init__(self)
        
        # check for explicit descriptor
        if kw.has_key('descr'):
            descr = kw.pop('descr')
            self._record.id = descr
            self._record.name = descr
            self._record.description = descr
        
        # define a set for uniq tags
        self._tags = set(self._record.annotations.setdefault('tags',[]))
    
    
    # define some simple interface to biopython internals
    
    @property
    def seq(self):
        return self._record.seq.tostring()
    
    @seq.setter
    def seq(self,s):
        self._record.seq = Seq(data=s,alphabet=NucleotideAlphabet())
    
    @property
    def descr(self):
        return self._record.descr
    
    @descr.setter
    def descr(self,d):
        self._record.id = d
        self._record.name = d
        self._record.description = d
    
    
    # define interface to tags object
    
    def add_tags(self,tags):
        if isinstance(tags,types.StringTypes):
            tags = [tags]
        elif isinstance(tags,types.ListType):
            tags = list(tags)
        else:
            raise TypeError, "value must be string type or list type"
        
        tags = set(tags)
        self._tags.update(tags)
        self._record.annotations['tags'] = list(self._tags)
        
        return self
    
    def add_tag(self,tag):
        return self.add_tags(tag)
    
    def has_tags(self,tags):
        if isinstance(tags,types.StringTypes):
            tags = [tags]
        elif isinstance(tags,types.ListType):
            tags = list(tags)
        else:
            raise TypeError, "value must be string type or list type"
        return set(tags) <= self._tags
    
    def has_tag(self,tag):
        return self.has_tags(tag)
    
    def del_tags(self,tags):
        if isinstance(tags,types.StringTypes):
            tags = [tags]
        elif isinstance(tags,types.ListType):
            tags = list(tags)
        else:
            raise TypeError, "value must be string type or list type"
        
        for tag in tags:
            self._tags.remove(tag)
        
        self._record.annotations['tags'] = list(self._tags)
        return self
    
    def del_tag(self,tag):
        return self.del_tags(tag)
    
    
    # define some functional interface:
    
    @property
    def junction(self):
        return self._record.features[self._features['CDR3-IMGT'][0]].extract(self._record.seq.tostring())
    
    @property
    def cdr3(self):
        return len(self.junction)
    
    @property
    def v(self):
        return self._record.features[self._features['V-REGION'][0]].qualifiers['allele']
    
    @property
    def d(self):
        return self._record.features[self._features['D-REGION'][0]].qualifiers['allele']
    
    @property
    def j(self):
        return self._record.features[self._features['J-REGION'][0]].qualifiers['allele']
    
    @property
    def vj(self):
        return '|'.join([self.v,self.j])
    
    @property
    def vdj(self):
        return '|'.join([self.v,self.d,self.j])
    
    def __len__(self):
        return len(self.seq)
    
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        return self._record.format('imgt')

# ================
# = INPUT/OUTPUT =
# ================

def parse_VDJXML(inputfile):
    """Parser for VDJXML
    
    Really just a wrapper around SeqIO
    """
    for record in SeqIO.parse(inputfile,'imgt'):
        yield ImmuneChain(record)

def filter_parse_VDJXML(inputfile,predicate):
    """Parser that takes a predicate function"""
    for record in SeqIO.parse(inputfile,'imgt'):
        chain = ImmuneChain(record)
        if predicate(chain):
            yield chain
