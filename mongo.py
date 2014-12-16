# Copyright 2014 Uri Laserson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pymongo

from seqtools import simplifySeqRecord, complicateSeqRecord
from vdj import ImmuneChain

def encode_chain(chain):
    document = simplifySeqRecord(chain)
    if hasattr(chain,'v'): document['v'] = chain.v
    if hasattr(chain,'d'): document['d'] = chain.d
    if hasattr(chain,'j'): document['j'] = chain.j
    if hasattr(chain,'junction_nt'): document['junction_nt'] = chain.junction_nt
    if hasattr(chain,'junction_aa'): document['junction_aa'] = chain.junction_aa
    if hasattr(chain,'num_mutations'): document['num_mutations'] = chain.num_mutations
    return document

def decode_document(document):
    assert document["__SeqRecord__"]
    return ImmuneChain(complicateSeqRecord(document))

def connect_to_spleen(connect_to='vaccination'):
    connection = pymongo.Connection('spleen.res.med.harvard.edu',27017)
    db = connection[connect_to]
    db.authenticate("mongodbuser","asdffdsa")
    return db

def connect_to_lymph(connect_to='vaccination'):
    connection = pymongo.Connection('hero1614.rc.fas.harvard.edu',27017)
    db = connection[connect_to]
    return db

def connect_to_localhost(connect_to='vaccination'):
    connection = pymongo.Connection('localhost',27017)
    db = connection[connect_to]
    return db


##############################################################################
# doesn't work if trying to insert an ImmuneChain only.  This will work if you
# try to insert a SON object that has an ImmuneChain as one of the values.  So
# you'd have to wrap the object in something like: {"chain": chain_obj}

from pymongo.son_manipulator import SONManipulator
from Bio.SeqRecord import SeqRecord

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
