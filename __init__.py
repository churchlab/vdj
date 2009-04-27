# __init__.py

# VDJ package: toolset for vdj manipulations

import numpy as np
from numpy import ma
import pylab
import scipy
import scipy.stats
import matplotlib
import seqtools
import refseq
import clones
import vdjexcept
import alignment
import analysis
from Bio.SeqRecord import SeqRecord
import cPickle
import operator
import countdata
from xml.sax.saxutils import escape, unescape
import xml.sax
import xml.sax.handler

spearmanr = scipy.stats.spearmanr
pearsonr  = scipy.stats.pearsonr

IGHV = refseq.IGHV
IGHD = refseq.IGHD
IGHJ = refseq.IGHJ

NV = len(IGHV)
ND = len(IGHD)
NJ = len(IGHJ)

#===============================================================================

# ===================
# = DATA STRUCTURES =
# ===================

class ImmuneChain(object):
	'''
		Data structure to represent an immune chain.
		
		VDJ is a triplet of strings representing the VDJ region that is used
		(computed elsewhere, however).
		
		Do not use commas in any internal descriptors
	'''
	
	def __init__(self, seq='', func='', v='', d='', j='', ighc='', cdr3=0, junction='', descr='', tags=set([])):
		'''
			seq is 5'->3'
			vdj is triplet of strings representing VDJ usage
			descr is optional additional identifier
		'''
		self.seq = seq.upper()
		self.descr = descr
		self.tags = set(tags)	# tag for sample number/experiment etc
		self.v = v
		self.d = d
		self.j = j
		self.ighc = ighc
		self.cdr3 = cdr3
		self.junction = junction
		self.func = func
	
	def __str__(self):
		return self.__repr__()
	
	def __repr__(self):
		return self.getXML()
	
	def getXML(self):
		xmlstring = ''
		xmlstring += '<ImmuneChain>\n'
		xmlstring += '\t<descr>' 	+ escape(self.descr) + 		'</descr>\n'
		xmlstring += '\t<seq>' 		+ escape(self.seq) + 		'</seq>\n'
		xmlstring += '\t<v>' 		+ escape(self.v) + 			'</v>\n' 
		xmlstring += '\t<d>' 		+ escape(self.d) + 			'</d>\n'
		xmlstring += '\t<j>' 		+ escape(self.j) + 			'</j>\n'
		xmlstring += '\t<ighc>' 	+ escape(self.ighc) + 		'</ighc>\n'
		xmlstring += '\t<cdr3>' 	+ escape(self.cdr3) + 		'</cdr3>\n'	# measured in nt
		xmlstring += '\t<junction>' + escape(self.junction) + 	'</junction>\n'
		xmlstring += '\t<func>' 	+ escape(self.func) + 		'</func>\n'
		for tag in self.tags:
			xmlstring += '\t<tag>' + escape(tag) + '</tag>\n'
		xmlstring += '</ImmuneChain>\n'
		return xmlstring


# ====================
# = Repertoire class =
# ====================

class Repertoire(object):
		
	def __init__( self, chains=[], metatags=[] ):
		'''
		must provide list of ImmuneChain objects
		they don't have to have full information
		'''
		if isinstance(chains,ImmuneChain):
			chains = [chains]
		
		if isinstance(chains,parseImmuneChains):
			chains = list(chains)
		
		# set global tags for repertoire
		self.metatags = set(metatags)
		
		# init primary datatype
		self.chains = np.array(chains,dtype=np.object)
		
		# collect all the tags in the set of ImmuneChains
		self.tags = {}
		
		# ensure that at least, all identifiers in refseq are present
		for ident in refseq.IGHV + refseq.IGHD + refseq.IGHJ:
			self.tags[ident] = []
		
		for (i,chain) in enumerate(self.chains):
			self.processTags(i,chain)
		
		self.uniqueifyTags()
		
		return
	
	
	# ======================
	# = Indexing/retrieval =
	# ======================
		
	def __getitem__(self,keys):
		'''
		get ImmuneChains out of this object
		takes anything that numpy arrays can take
		'''
		# num dim
		if not isinstance(keys,list) and not isinstance(keys,tuple):
			keys = (keys,)
		
		if isinstance(keys[0],int) or isinstance(keys[0],slice):
			keys = list(keys)
			if len(keys)==1 and isinstance(keys[0],int):
				return self.chains[keys[0]]
			if len(keys)==1 and isinstance(keys[0],slice):
				return Repertoire(self.chains[keys[0]],self.metatags)
			
			return Repertoire(self.chains[keys],self.metatags)
		
		raise IndexError, "you must've indexed incorrectly because you shouldn't be here"
		
	def get_chains_AND(self,args):
		idxs = self.get_idxs_AND(args)
		return Repertoire(self.chains[idxs],self.metatags)
	
	def get_chains_OR(self,args):
		idxs = self.get_idxs_OR(args)		
		return Repertoire(self.chains[idxs],self.metatags)
		
	def get_idxs_AND(self,args):
		identifiers = self.tags.keys()	  # all valid identifiers
		for key in args:	# check for all valid identifiers
			if key not in identifiers:
				print 'WARNING: ' + key + ' is not a recognized identifier'
				return []
		
		idxs = range(len(self.chains))
		for key in args:
			idxs = list(set(idxs) & set(self.tags[key]))
		
		return idxs
	
	def get_idxs_OR(self,args):
		args = list(args)	# necessary bc of pop operation below
		identifiers = self.tags.keys()
		for (i,key) in enumerate(args):
			if key not in identifiers:
				print 'WARNING: ' + key + ' is not a recognized identifier'
				args.pop(i)
				
		idxs = set([])
		for key in args:
			idxs.update(self.tags[key])
		idxs = list(idxs)
		
		return idxs
	
	def get_idxs_fullVDJ(self):
		# indexing of reference lists starts at 1 bc the first elt is ''
		Vs = set(self.get_idxs_OR(refseq.IGHV[1:]))
		Ds = set(self.get_idxs_OR(refseq.IGHD[1:]))
		Js = set(self.get_idxs_OR(refseq.IGHJ[1:]))
		idxs = list(Vs & Ds & Js)
		idxs.sort()
		return idxs
	
	def get_idxs_fullVJ(self):
		Vs = set(self.get_idxs_OR(refseq.IGHV[1:]))
		Js = set(self.get_idxs_OR(refseq.IGHJ[1:]))
		idxs = list(Vs & Js)
		idxs.sort()
		return idxs
	
	def get_idxs_fullVJCDR3(self):
		VJs   = set(self.get_idxs_fullVJ())
		CDR3s = set([i for (i,chain) in enumerate(self.chains) if chain.cdr3 % 3 == 0 and chain.cdr3 > 0])
		idxs = list(VJs & CDR3s)
		idxs.sort()
		return idxs		
	
	def get_chains_fullVDJ(self):
		return self[ self.get_idxs_fullVDJ() ]
	
	def get_chains_fullVJ(self):
		return self[ self.get_idxs_fullVJ() ]
	
	def get_chains_fullVJCDR3(self):
		return self[ self.get_idxs_fullVJCDR3() ]
	
	# ==================================
	# = Appending/extending repertoire =
	# ==================================
	
	def __add__(self,other):
			'''
			combine two repertoires
			'''
			# chain objects must be concatenated
			# tags object must be updated
			return Repertoire(np.append(self.chains,other.chains), (self.metatags | other.metatags) )
	
	def __radd__(self,other):
		return Repertoire.__add__(other,self)
	
	def __iadd__(self,other):
		for chain in other:
			self.append(chain)
		return self
	
	def __len__(self):
		return len(self.chains)
	
	def append(self,chain):
		i = len(self.chains)
		self.chains = np.append(self.chains,chain)
		self.processTags(i,chain)
		
		# not necessary to uniqueify the tags as each
		# time i add another one it's a higher unrepresented
		# number
		
		return self
		
	def extend(self,rep):
		return self.__iadd__(rep)
	
	# ======================
	# = Iterator interface =
	# ======================
	
	def __iter__(self):
		return self.chains.__iter__()
	
	# =============
	# = Utilities =
	# =============
	
	def add_tags(self,tagset):
		if isinstance(tagset,str): tagset = [tagset]
		for (i,chain) in enumerate(self.chains):
			chain.tags.update(tagset)
			for tag in tagset:
				try: self.tags[tag] += [i]
				except KeyError,e: self.tags[tag] = [i]
		return
	
	def del_tags(self,tagset):
		if isinstance(tagset,str): tagset = [tagset]
		for tag in tagset:
			if tag not in self.tags.keys():
				print 'WARNING: ' + tag + 'is not a recognized tag'
				continue
			idxs = self.tags[tag]
			for chain in self.chains[idxs]
				chain.tags.remove(tag)
			del self.tags[tag]
		return
	
	def add_metatags(self,metatagset):
		if isinstance(metatagset,str): metatagset = [metatagset]
		self.metatags.update(metatagset)
		return
	
	def del_metatags(self,metatagset):
		if isinstance(metatagset,str): metatagset = [metatagset]
		for tag in metatagset:
			if tag not in self.metatags:
				print 'WARNING: ' + tag + 'is not a recognized tag'
				continue
			self.metatags.remove(tag)
		return
		
	def processTags(self,i,chain):
		# given ImmuneChain and position in array, get tags and add
		# them to the tags dict
		for tag in chain.tags:
			try: self.tags[tag] += [i]
			except KeyError,e: self.tags[tag] = [i]
		
		try: self.tags[chain.v] += [i]
		except KeyError,e: self.tags[chain.v] = [i]
		
		try: self.tags[chain.d] += [i]
		except KeyError,e: self.tags[chain.d] = [i]
		
		try: self.tags[chain.j] += [i]
		except KeyError,e: self.tags[chain.j] = [i]
		
		try: self.tags[chain.descr] += [i]
		except KeyError,e: self.tags[chain.descr] = [i]
		
		return
	
	def uniqueifyTags(self):
		for tag in self.tags.keys():
			self.tags[tag] = list(set(self.tags[tag]))
		return
	
	def __str__(self):
		return self.__repr__()
	
	def __repr__(self):
		return self.getXML()
	
	def getXML(self):
		xmlstring = ''
		xmlstring += '<Repertoire>\n\n'
		xmlstring += '<Meta>\n'
		for tag in self.metatags:
			xmlstring += '\t<metatag>' + escape(tag) + '</metatag>\n'
		xmlstring += '</Meta>\n\n'
		xmlstring += '<Data>\n'
		for chain in self.chains:
			xmlstring += chain.getXML() + '\n'
		xmlstring += '</Data>\n\n'
		xmlstring += '</Repertoire>\n'
		return xmlstring

#===============================================================================

# ================
# = Input/Output =
# ================

class VDJXMLhandler(xml.sax.handler.ContentHandler):
	def __init__(self,dataobject):
		# container for data
		self.data = dataobject
		if isinstance(self.data,Repertoire):
			self.mode = 'Repertoire'
		elif isinstance(self.data,list):
			self.mode = 'ImmuneChain'
		else:
			raise TypeError, 'you provided the wrong type of container to VDJXMLhandler'
		
		self.buffer = ''
		self.numChains = 0
	
	def startElement(self,name,attrs):
		elif name == 'ImmuneChain':
			self.currentChain = ImmuneChain()
	
	def characters(self,content):
		if saveData == True:
			self.buffer += unescape(content)
	
	def endElement(self,name):
		if name == 'metatag':
			if mode == 'Repertoire':
				self.data.add_metatags(self.buffer)
		elif name == 'ImmuneChain':
			self.data.append(self.currentChain)
			self.numChains += 1
		elif name == 'descr':
			self.currentChain.descr = self.buffer
		elif name == 'seq':
			self.currentChain.seq = self.buffer
		elif name == 'v':
			self.currentChain.v = self.buffer
		elif name == 'd':
			self.currentChain.d = self.buffer
		elif name == 'j':
			self.currentChain.j = self.buffer
		elif name == 'ighc':
			self.currentChain.ighc = self.buffer
		elif name == 'cdr3':
			self.currentChain.cdr3 = eval(self.buffer)
		elif name == 'junction':
			self.currentChain.junction = self.buffer
		elif name == 'func':
			self.currentChain.func = self.buffer
		elif name == 'tag':
			self.currentChain.tags.update([self.buffer])
		else:
			self.buffer = ''
		self.buffer = ''
	
	def endDocument(self):
		print "Mode is: " + self.mode
		if mode == 'Repertoire':
			print "Repertoire Metatags:"
			for tag in self.data.metatags:
				print tag
		print "Number of ImmuneChain objects read: " + str(self.numChains)

def readVDJ(inputfile,mode='Repertoire'):
	'''
	function to load a Repertoire and return it in a Repertoire object
	or a list of ImmuneChains
	'''
	if isinstance(inputfile,str):
		ip = open(inputfile,'r')
	elif isinstance(inputfile,file):
		ip = inputfile
	
	if mode == 'Repertoire':
		data = Repertoire()
	elif mode == 'ImmuneChain':
		data = []
	handler = VDJXMLhandler(data)
	saxparser = xml.sax.make_parser()
	saxparser.setContentHandler(handler)
	saxparser.parse(ip)
	
	if isinstance(inputfile,str):
		ip.close()
	
	return data

def writeVDJ(data, handle):
	print >>handle, '<?xml version="1.0"?>'
	if isinstance(data,Repertoire):
		print >>handle, data
	elif isinstance(data,list) and isinstance(data[0],ImmuneChain):
		for chain in data:
			print >>handle, chain

#===============================================================================

# ============
# = Counting =
# ============

def countsVJ(rep):
	cn = np.zeros( (NV,NJ) )
	for chain in rep.chains:
		cn[refseq.IGHVdict[chain.v],refseq.IGHJdict[chain.j]] += 1
	return cn

def countsVJ_1D(rep):
	return countsVJ(rep).ravel()

def countsVDJ(rep):
	cn = np.zeros( (NV,ND,NJ) )
	for chain in rep.chains:
		cn[refseq.IGHVdict[chain.v],refseq.IGHDdict[chain.d],refseq.IGHJdict[chain.j]] += 1
	return cn

def countsVDJ_2D(rep):
	cn = countsVDJ(rep)
	return cn.reshape(NV,ND*NJ)

def countsVDJ_1D(rep):
	return countsVDJ(rep).ravel()

def countsVJCDR3(rep,cdrlow=3,cdrhigh=99):
	numlens = (cdrhigh-cdrlow)/3 + 1
	cn = np.zeros( (NV,NJ,numlens) )
	for chain in rep.chains:
		if chain.cdr3 >= cdrlow and chain.cdr3 <= cdrhigh and chain.cdr3 % 3 == 0:
			cn[refseq.IGHVdict[chain.v],refseq.IGHJdict[chain.j],(chain.cdr3-cdrlow)/3] += 1
	return cn

def countsVJCDR3_2D(rep,cdrlow=3,cdrhigh=99):
	cn = countsVJCDR3(rep,cdrlow,cdrhigh)
	return cn.reshape(NV,cn.shape[1]*cn.shape[2])

def countsVJCDR3_1D(rep,cdrlow=3,cdrhigh=99):
	return countsVJCDR3(rep,cdrlow,cdrhigh).ravel()
