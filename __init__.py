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

spearmanr = scipy.stats.spearmanr
pearsonr  = scipy.stats.pearsonr

IGHV = refseq.IGHV
IGHD = refseq.IGHD
IGHJ = refseq.IGHJ

NV = len(IGHV)
ND = len(IGHD)
NJ = len(IGHJ)

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
		self.tags = tags	# tag for sample number/experiment etc
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
		callstring = 'ImmuneChain{'
		callstring += self.seq + ','
		callstring += self.func + ','
		callstring += self.v + ',' + self.d + ',' + self.j + ','
		callstring += self.ighc + ','
		callstring += str(self.cdr3) + ','
		callstring += self.junction + ','
		callstring += self.descr
		for tag in self.tags: callstring += ',' + tag
		callstring += '}'
		return callstring

class parseImmuneChains(object):
	'''
	returns iterator for use in parsing output of different aligners
	'''
	
	def __init__(self,handle,format='ImmuneChain'):
		self.ip = handle
		if format == 'ImmuneChain': self.next_impl = self.read_immunechain
		elif format == 'abacus': self.next_impl = self.read_abacus		# less important
		elif format == 'V-QUEST': self.next_impl = self.read_v_quest	# less important
	
	def __iter__(self):
		return self
	
	def next(self):
		return self.next_impl()
	
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


# ====================
# = Repertoire class =
# ====================

class Repertoire(object):
		
	def __init__( self, chains=[] ):
		'''
		must provide list of ImmuneChain objects
		they don't have to have full information
		'''
		if isinstance(chains,ImmuneChain):
			chains = [chains]
		
		if isinstance(chains,parseImmuneChains):
			chains = list(chains)
		
		# init primary datatype
		self.chains = np.array(chains,dtype=np.object)
		
		# collect all the tags in the set of ImmuneChains
		tags = {}
		
		# ensure that at least, all identifiers in refseq are present
		for ident in refseq.IGHV + refseq.IGHD + refseq.IGHJ:
			tags[ident] = []
		
		for (i,chain) in enumerate(chains):
			for tag in chain.tags:
				try: tags[tag] += [i]
				except KeyError,e: tags[tag] = [i]
			
			try: tags[chain.v] += [i]
			except KeyError,e: tags[chain.v] = [i]
			
			try: tags[chain.d] += [i]
			except KeyError,e: tags[chain.d] = [i]
			
			try: tags[chain.j] += [i]
			except KeyError,e: tags[chain.j] = [i]
			
			try: tags[chain.descr] += [i]
			except KeyError,e: tags[chain.descr] = [i]
		
		# unique-ify the lists for all tags
		for tag in tags.keys():
			tags[tag] = list(set(tags[tag]))
		
		self.tags = tags
		
		return
	
	# ======================
	# = Indexing/retrieval =
	# ======================
		
	def __getitem__(self,keys):
		'''
		get ImmuneChains out of this object
		
		if given a tuple of strings, it returns all chains that
		match ALL identifiers in the tuple
		the string can also match the tags that identify ImmuneChains
		'''
		# num dim
		if not isinstance(keys,list) and not isinstance(keys,tuple):
			keys = (keys,)
		
		# convert to string interface if necessary
		if isinstance(keys[0],int) or isinstance(keys[0],slice):
			keys = list(keys)
			
			# if obj[:]
			if len(keys)==1 and isinstance(keys[0],int):
				return self.chains[keys[0]]
			if len(keys)==1 and isinstance(keys[0],slice):
				return Repertoire(self.chains[keys[0]])
			
			return Repertoire(self.chains[keys])
			
		if isinstance(keys[0],str):
			# process keys
			identifiers = self.tags.keys()	  # all valid identifiers
			for key in keys:	# check for all valid identifiers
				if key not in identifiers:
					raise KeyError, key+' is not a recognized identifier'
			
			idxs = range(len(self.chains))
			for key in keys:
				idxs = list(set(idxs) & set(self.tags[key]))
			
			return Repertoire(self.chains[idxs])
		
		raise IndexError, "you must've indexed incorrectly because you shouldn't be here"
		
	def get_chains_AND(self,args):
		return self.__getitem__(args)
	
	def get_chains_OR(self,args):
		identifiers = self.tags.keys()
		for key in args:
			if key not in identifiers:
				raise KeyError, key+' is not a recognized identifier'
		
		idxs = set([])
		for key in args:
			idxs.update(self.tags[key])
		idxs = list(idxs)
		
		return Repertoire(self.chains[idxs])
		
	def get_idxs_AND(self,args):
		identifiers = self.tags.keys()	  # all valid identifiers
		for key in args:	# check for all valid identifiers
			if key not in identifiers:
				raise KeyError, key+' is not a recognized identifier'
		
		idxs = range(len(self.chains))
		for key in args:
			idxs = list(set(idxs) & set(self.tags[key]))
		
		return idxs
	
	def get_idxs_OR(self,args):
		identifiers = self.tags.keys()
		for key in args:
			if key not in identifiers:
				raise KeyError, key+' is not a recognized identifier'
		
		idxs = set([])
		for key in args:
			idxs.update(self.tags[key])
		idxs = list(idxs)
		
		return idxs
	
	def get_idxs_fullVDJ(self):
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
	
	def get_rep_fullVDJ(self):
		return self[ self.get_idxs_fullVDJ() ]
	
	def get_rep_fullVJ(self):
		return self[ self.get_idxs_fullVJ() ]
	
	def get_rep_fullVJCDR3(self):
		return Repertoire([chain for chain in self if chain.v in refseq.IGHV[1:] and chain.j in refseq.IGHJ[1:] and chain.cdr3 % 3 == 0 and chain.cdr3 > 0])
	
	# ==================================
	# = Appending/extending repertoire =
	# ==================================
	
	def __add__(self,other):
			'''
			combine two repertoires
			'''
			# chain objects must be concatenated
			# tags object must be updated
			return Repertoire(np.append(self.chains,other.chains))
	
	def __radd__(self,other):
		return Repertoire.__add__(other,self)
	
	def __iadd__(self,other):
		n = len(self.chains)
		m = len(other.chains)
		self.chains = np.append(self.chains,other.chains)
		
		tags = self.tags
		for (j,chain) in enumerate(other.chains):
			i = n+j
			for tag in chain.tags:
				try: tags[tag] += [i]
				except KeyError,e: tags[tag] = [i]
			
			try: tags[chain.v] += [i]
			except KeyError,e: tags[chain.v] = [i]
			
			try: tags[chain.d] += [i]
			except KeyError,e: tags[chain.d] = [i]
			
			try: tags[chain.j] += [i]
			except KeyError,e: tags[chain.j] = [i]
			
			try: tags[chain.descr] += [i]
			except KeyError,e: tags[chain.descr] = [i]
		
		# unique-ify the lists for all tags
		for tag in tags.keys():
			tags[tag] = list(set(tags[tag]))
		
		self.tags = tags
		
		return self
	
	def __len__(self):
		return len(self.chains)
	
	def append(self,chain):
		return self.__iadd__(self,Repertoire(chain))
	
	def extend(self,rep):
		return self.__iadd__(rep)
	
	# ======================
	# = Iterator interface =
	# ======================
	
	def __iter__(self):
		return self.chains.__iter__()
	
	# ============
	# = Counting =
	# ============
	
	def countsVJ(self):
		cn = np.zeros( (NV,NJ) )
		for chain in self.chains:
			cn[refseq.IGHVdict[chain.v],refseq.IGHJdict[chain.j]] += 1
		return cn
	
	def countsVJ_1D(self):
		return self.countsVJ().ravel()
	
	def countsVDJ(self):
		cn = np.zeros( (NV,ND,NJ) )
		for chain in self.chains:
			cn[refseq.IGHVdict[chain.v],refseq.IGHDdict[chain.d],refseq.IGHJdict[chain.j]] += 1
		return cn
	
	def countsVDJ_2D(self):
		cn = self.countsVDJ()
		return cn.reshape(NV,ND*NJ)
	
	def countsVDJ_1D(self):
		return self.countsVDJ().ravel()
	
	def countsVJCDR3(self,cdrlow=3,cdrhigh=99):
		numlens = (cdrhigh-cdrlow)/3 + 1
		cn = np.zeros( (NV,NJ,numlens) )
		for chain in self.chains:
			if chain.cdr3 >= cdrlow and chain.cdr3 <= cdrhigh and chain.cdr3 % 3 == 0:
				cn[refseq.IGHVdict[chain.v],refseq.IGHJdict[chain.j],(chain.cdr3-cdrlow)/3] += 1
		return cn
	
	def countsVJCDR3_2D(self,cdrlow=3,cdrhigh=99):
		cn = self.countsVJCDR3(cdrlow,cdrhigh)
		return cn.reshape(NV,cn.shape[1]*cn.shape[2])
	
	def countsVJCDR3_1D(self,cdrlow=3,cdrhigh=99):
		return self.countsVJCDR3(cdrlow,cdrhigh).ravel()
	
	# =============
	# = Utilities =
	# =============
	
	def set_tags(self,tagset):
		for (i,chain) in enumerate(self.chains):
			chain.tags.update(tagset)
			for tag in tagset:
				try: self.tags[tag] += [i]
				except KeyError,e: self.tags[tag] = [i]
		return
	
	def __str__(self):
		return self.__repr__()
	
	def __repr__(self):
		strrep = ''
		for chain in self.chains:
			strrep += str(chain) + '\n'
		return strrep
	
