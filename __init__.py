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

class ImmuneChain(object):
    """Data structure to represent an immune chain.
    
    The underlying data structure is actually a biopython `SeqRecord` object.
    This class wraps it in a way that maintains the simple-to-use interface to
    get at some common annotations. It also knows how to print out it's data
    as IMGT-flavored INSDC (GenBank/EMBL).
    
    Sequences must be nucleotides, not protein.
    """
    
    def __init__(self,**kw):
        """Initialize ImmuneChain
        
        seq is 5'->3'
        """
        # first we define our underlying SeqRecord object
        if 'record' in kw:  # init with SeqRecord object
            self._record = kw['record']
        else:   # otherwise, initialize manually
            if 'seq' in kw:
                seq = Seq(data=kw['seq'],alphabet=NucleotideAlphabet())
                self._record = SeqRecord(seq=seq,id='',name='',description='')
            else:
                self._record = SeqRecord(seq=UnknownSeq(0,alphabet=NucleotideAlphabet()),id='',name='',description='')
        
        # check for explicit descriptor
        if kw.has_key('descr'):
            descr = kw.pop('descr')
            self._record.id = descr
            self._record.name = descr
            self._record.description = descr
        
        # define dictionary of features for faster lookup
        self._features = {}
        for (i,feature) in enumerate(self_record.features):
            self._features.setdefault(feature.type,[]).append(i)
        
        # define a set for uniq tags
        self._tags = set(self._record.annotations.setdefault('tags',[]))
        
        # define list of special attribute names. These attributes are all
        # handled wrt the underlying SeqRecord object object automatically. if i
        # try to set/get any other attribute, i must catch it to ensure its info
        # gets encoded in the SeqRecord. this is done in the __setattr__ and
        # __getattr__ fns
        self._reserved = [ 'seq',
                           'descr',
                           'tags',
                           'cdr3',
                           'junction',
                           'v',
                           'd',
                           'j',
                           'vj',
                           'vdj' ]
    
    # define some simple interface to biopython internals
    
    def __getattr__(self,name):
        # This function should only get called if I am looking for an attribute that
        # didn't already have a setter defined or a default method.  In this case, I
        # search the annotations dictionary of the underlying SeqRecord to try to find
        # the information.
        
        # if name not in self._reserved:
        return self._record.annotations[name]
        # return object.__getattr__(self,name)
    
    def __setattr__(self,name,value):
        object.__setattr__(self,name,value)         # python object attr setting
        if name not in self._reserved:
            self._record.annotations[name] = value  # reflect attr in SeqRecord annotations
    
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
        self._record.descr = d
    
    def set_annot(self,name,value):
        if name in self._reserved:
            raise ValueError: "'%s' annotation is reserved" % name
        
        if isinstance(value,types.StringTypes):
            value = str(value)
        elif isinstance(value,types.ListType):
            value = list(value)
        else:
            raise TypeError, "value must be string type or list type"
        
        self._record.annotations[name] = value
        return self
    
    def has_annot(self,name):
        return name in self._record.annotations
    
    def del_annot(self,name):
        if name == 'tags':
            raise ValueError: "'tags' annotation is reserved"
        
        if self.has_annot(name):
            self._record.annotations.pop(name)
        else:
            raise KeyError, ("%s is not an annotation" % name)
        return self
    
    def add_feature(self,start=None,end=None,type='',strand=None,qualifiers=None):
        if start == None and end == None:
            raise ValueError, "if there is no location, use an annotation"
        
        if start != None and end != None:
            location = FeatureLocation(start,end)
        elif start != None and end == None:
            location = FeatureLocation(start,start)
        elif start == None and end != None:
            location = FeatureLocation(end,end)
        
        feature = SeqFeature(location=location,type=type,strand=strand,qualifiers=qualifiers)
        
        self._record.features.append(feature)
        self._features.setdefault(feature.type,[]).append(len(self._record.features) - 1)
        return self
    
    def has_feature(self,type):
        return type in self._features
    
    def del_feature(self,type):
        idxs = self._features.pop(type)
        idxs.sort(reverse=True)
        for i in idxs:
            self._record.features.pop(i)
        return self
    
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
    
    def has_tags(self,tags):
        if isinstance(tags,types.StringTypes):
            tags = [tags]
        elif isinstance(tags,types.ListType):
            tags = list(tags)
        else:
            raise TypeError, "value must be string type or list type"
        return set(tags) <= self._tags
    
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
    
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        return self._record.format('imgt')
    
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
    
    @property
    def tags(self):
        return self._tags
    
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
