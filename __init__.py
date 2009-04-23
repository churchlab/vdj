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
		for (i,chain) in enumerate(self.chains):
			chain.tags.update(tagset)
			for tag in tagset:
				try: self.tags[tag] += [i]
				except KeyError,e: self.tags[tag] = [i]
		return
	
	def del_tags(self,tagset):
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
		self.metatags.update(metatagset)
		return
	
	def del_metatags(self,metatagset):
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



class parseImmuneChains(object):
	'''
	returns iterator for use in parsing output of different aligners
	'''
	
	def __init__(self,handle,format='VDJXML'):
		self.ip = handle
		if format   == 'VDJXML':
			self.next_impl = self.read_repXML
			handler = RepertoireXMLhandler()
			saxparser = xml.sax.make_parser()
			saxparser.setContentHandler(handler)
		elif format == 'ImmuneChain': self.next_impl = self.read_immunechain
		elif format == 'abacus': self.next_impl = self.read_abacus		# less important
		elif format == 'V-QUEST': self.next_impl = self.read_v_quest	# less important
	
	def __iter__(self):
		return self
	
	def next(self):
		return self.next_impl()
	
	
	class RepertoireXMLhandler(xml.sax.handler.ContentHandler):
		def startElement(self,name,attrs):
			if name == 'metatag':
				self.inMetatag = True
				self.saveData = True
			elif name == 'ImmuneChain':
				self.inImmuneChain = True
				self.saveData = False
				
		def characters(self,content):
			pass

		def endElement(self,name):
			if name == 'metatag':
				self.inMetatag = False
				self.saveData = False
	
	def read_repXML(self):
		
	
	def read_immunechain(self):
		# clean data
		line = self.ip.next().strip()
		data = line.strip('ImmuneChain').strip('{}').split(',')
		numtags = len(data) - 9		# there are 9 req'd fields before tags are included
		
		# build ImmuneChain
		return ImmuneChain( seq=data[0],\
							func=data[1],\
							v=data[2],\
							d=data[3],\
							j=data[4],\
							ighc=data[5],\
							cdr3=eval(data[6]),\
							junction=data[7],\
							descr=data[8],\
							tags=set(data[9:]) )
	
	def read_abacus(self):
		line = self.ip.next()
		name = line.split('|')[0].strip()
		rawvdj = line.split('|')[1:]
		v = ''
		d = ''
		j = ''
		for segment in rawvdj:
			if segment.split('*')[0][3] == 'V':
				v = segment.split('*')[0]
			elif segment.split('*')[0][3] == 'D':
				d = segment.split('*')[0]
			elif segment.split('*')[0][3] == 'J':
				j = segment.split('*')[0]
		parsedvdj = ( v, d, j )
		return ImmuneChain(descr = name, vdj = parsedvdj)
	
	def read_v_quest(self):
		line = self.ip.next()
		rawdata = line.split(',')
		return ImmuneChain( seq		 = rawdata[0][13:], \
							func	 = rawdata[1], \
							vdj		 = (rawdata[2][2:-1],rawdata[3][2:-1],rawdata[4][2:-2]), \
							cdr3	 = eval(rawdata[5]), \
							junction = rawdata[6], \
							descr	 = rawdata[7][:-4] )	# for the ');\n'

def writeImmuneChains(chains, handle):
	for chain in chains:
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