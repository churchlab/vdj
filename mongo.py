from pymongo.son_manipulator import SONManipulator
from Bio.SeqRecord import SeqRecord

from seqtools import simplifySeqRecord, complicateSeqRecord
from vdj import ImmuneChain

encode_chain = simplifySeqRecord
def decode_chain(document):
    assert document["__SeqRecord__"]
    return ImmuneChain(complicateSeqRecord(document))






##############################################################################
# doesn't work if trying to insert an ImmuneChain only.  This will work if you
# try to insert a SON object that has an ImmuneChain as one of the values.  So
# you'd have to wrap the object in something like: {"chain": chain_obj}

class ImmuneChainTransform(SONManipulator):
    def transform_incoming(self,son,collection):
        for (key,value) in son.items():
            if isinstance(value,SeqRecord):
                son[key] = simplifySeqRecord(value)
            elif isinstance(value,dict):
                son[key] = self.transform_incoming(value,collection)
        return son
    
    def transform_outgoing(self,son,collection):
        for (key,value) in son.items():
            if isinstance(value,dict):
                if "__SeqRecord__" in value and value["__SeqRecord__"]:
                    son[key] = ImmuneChain(complicateSeqRecord(value))
                else:
                    son[key] = self.transform_outgoing(value,collection)
        return son
