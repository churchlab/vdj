import types

from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import NucleotideAlphabet
from Bio import SeqIO

# select which IMGT reference set to load; human by default
organism = 'human'
# organism = 'mouse'

# ===================
# = DATA STRUCTURES =
# ===================

class ImmuneChain(SeqRecord):
    """Data structure to represent an immune chain.
    
    It extends a biopython SeqRecord object with some simpler interface for
    common analyses.
    """
    
    def __init__(self, *args, **kw):
        """Initialize ImmuneChain
        
        This is performed either with a prebuilt SeqRecord object or as a
        native SeqRecord object.
        """
        if isinstance(args[0],SeqRecord):   # pre-build SeqRecord
            self.init_with_SeqRecord(args[0])
        elif kw.has_key('record'):          # pre-built SeqRecord
            self.init_with_SeqRecord(kw['record'])
        else:   # native SeqRecord init
            SeqRecord.__init__(self,*args,**kw)
        
        # define a set for uniq tags
        self._tags = set(self.annotations.setdefault('tags',[]))
        
        # precompute hash on features for performance
        self.update_feature_dict()
    
    def init_with_SeqRecord(self,record):
        # Initialize self using existing SeqRecord object
        SeqRecord.__init__(self, seq=record.seq, id=record.id,
                            name=record.name, description=record.description,
                            dbxrefs=record.dbxrefs, features=record.features,
                            annotations=record.annotations,
                            letter_annotations=record.letter_annotations)
    
    def update_feature_dict(self):
        self._features = {}
        for (i,feature) in enumerate(self.features):
            self._features.setdefault(feature.type,[]).append(i)
    
    
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
        self.annotations['tags'] = list(self._tags)
        
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
        
        self.annotations['tags'] = list(self._tags)
        return self
    
    def del_tag(self,tag):
        return self.del_tags(tag)
    
    
    # define some functional interface:
    
    @property
    def junction(self):
        return self.features[self._features['CDR3-IMGT'][0]].extract(self._record.seq.tostring())
    
    @property
    def cdr3(self):
        return len(self.junction)
    
    @property
    def v(self):
        return self.features[self._features['V-REGION'][0]].qualifiers['allele']
    
    @property
    def d(self):
        return self.features[self._features['D-REGION'][0]].qualifiers['allele']
    
    @property
    def j(self):
        return self.features[self._features['J-REGION'][0]].qualifiers['allele']
    
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
        return self.format('imgt')

# ================
# = INPUT/OUTPUT =
# ================

def parse_imgt(inputfile):
    """Parser for VDJXML
    
    Really just a wrapper around SeqIO
    """
    for record in SeqIO.parse(inputfile,'imgt'):
        yield ImmuneChain(record)

def filter_parse_imgt(inputfile,predicate):
    """Parser that takes a predicate function"""
    for record in SeqIO.parse(inputfile,'imgt'):
        chain = ImmuneChain(record)
        if predicate(chain):
            yield chain
