import types
import xml.sax
import xml.sax.handler
import time
import operator

import numpy as np

import seqtools
import refseq
import alignmentcore



# ===================
# = DATA STRUCTURES =
# ===================

class ImmuneChain(object):
	"""Data structure to represent an immune chain."""
	
	def __init__(self, seq='', func='', v='', d='', j='', ighc='', cdr3=0, junction='', descr='', tags=set([])):
		"""Initialize ImmuneChain
		
		seq is 5'->3'
		"""
		self.seq = seq.upper()
		self.descr = descr
		if isinstance(tags,types.StringTypes): tags = [tags]
		self.tags = set(tags)	# tag for sample number/experiment etc
		self.v = v
		self.d = d
		self.j = j
		self.ighc = ighc
		self.cdr3 = cdr3
		self.junction = junction
		self.func = func
	
	# REWORK CDR3 AS PROPERTY LATER
	# def getCDR3len(self):
	# 	return len(self.junction)
	# 
	# def setCDR3len(self):
	# 	# empty method.  checks if set value corresponds 
	# 
	# cdr3 = property(getCDR3len)
	
	def add_tags(self,tagset):
		if isinstance(tagset,types.StringTypes): tagset = [tagset]
		self.tags.update(tagset)
	
	def __len__(self):
		return len(self.seq)
	
	def __str__(self):
		return self.__repr__()
	
	def __repr__(self):
		return self.getXML()
	
	def getXML(self):
		xmlstring = ''
		xmlstring += '<ImmuneChain>\n'
		xmlstring += '\t<descr>' 	+ self.descr + 		'</descr>\n'
		xmlstring += '\t<seq>' 		+ self.seq + 		'</seq>\n'
		xmlstring += '\t<v>' 		+ self.v + 			'</v>\n' 
		xmlstring += '\t<d>' 		+ self.d + 			'</d>\n'
		xmlstring += '\t<j>' 		+ self.j + 			'</j>\n'
		xmlstring += '\t<ighc>' 	+ self.ighc + 		'</ighc>\n'
		xmlstring += '\t<cdr3>' 	+ str(self.cdr3) + 	'</cdr3>\n'	# measured in nt
		xmlstring += '\t<junction>' + self.junction + 	'</junction>\n'
		xmlstring += '\t<func>' 	+ self.func + 		'</func>\n'
		for tag in self.tags:
			xmlstring += '\t<tag>' + tag + '</tag>\n'
		xmlstring += '</ImmuneChain>\n'
		return xmlstring

# ====================
# = Repertoire class =
# ====================

class Repertoire(object):
	"""Data structure to contain and access many ImmuneChain objects"""
		
	def __init__( self, chains=[], metatags=[] ):
		"""Initialize Repertoire object.  Similar to list, but VDJ enhanced functionality.
		
		chains -- a list of ImmuneChain objects (they don't require full info)
		
		underlying data structure is numpy list
		
		"""
		if isinstance(chains,ImmuneChain):
			chains = [chains]
		
		# set global tags for repertoire
		if isinstance(metatags,types.StringTypes): metatags = [metatags]
		self.metatags = set(metatags)
		
		# init primary datatype
		self.chains = np.array(chains,dtype=np.object)
		
		# collect all the tags in the set of ImmuneChains
		self.tags = {}
		
		# ensure that at least, all identifiers in refseq are present
		for ident in refseq.ALL_IDs[1:]:
			self.tags[ident] = []
		
		for (i,chain) in enumerate(self.chains):
			self.processTags(i,chain)
		
		self.uniqueifyTags()
		
		return
	
	
	# ======================
	# = Indexing/retrieval =
	# ======================
		
	def __getitem__(self,keys):
		"""Return subset of Repertoire using Numpy-like indexing
		
		If keys is just an integer, then return the ImmuneChain object.
		If keys is anything else, it will return a Repertoire object with
			the corresponding chains in it.
			
		The metatags are inherited too.
		
		"""
		
		# num dim
		if not isinstance(keys,list) and not isinstance(keys,tuple) and not isinstance(keys,np.ndarray):
			keys = (keys,)
		
		if len(keys) == 0:
			return Repertoire([],self.metatags)
		
		if isinstance(keys[0],int) or isinstance(keys[0],slice):
			keys = list(keys)
			if len(keys)==1 and isinstance(keys[0],int):
				return self.chains[keys[0]]
			if len(keys)==1 and isinstance(keys[0],slice):
				return Repertoire(self.chains[keys[0]],self.metatags)
			
			return Repertoire(self.chains[keys],self.metatags)
		
		raise IndexError, "you must've indexed incorrectly because you shouldn't be here"
		
	def get_chains_AND(self,args):
		"""Return subRepertoire that matches ALL strings in args"""
		return self[ self.get_idxs_AND(args) ]
	
	def get_chains_OR(self,args):
		"""Return subRepertoire that matches ANY strings in args"""
		return self[ self.get_idxs_OR(args) ]
		
	def get_idxs_AND(self,args):
		"""Return indexes of chains that match ALL strings in args"""
		if isinstance(args,types.StringTypes): args = [args]
		
		# check that requested identifiers are present
		identifiers = self.tags.keys()	  # all valid identifiers
		for key in args:	# check for all valid identifiers
			if key not in identifiers:
				print 'WARNING: ' + key + ' is not a recognized identifier.  Returning empty list.'
				return []
		
		idxs = range(len(self.chains))
		for key in args:
			idxs = list(set(idxs) & set(self.tags[key]))
		
		return idxs
	
	def get_idxs_OR(self,args):
		"""Return indexes of chains that match ANY strings in args"""
		if isinstance(args,types.StringTypes): args = [args]
		args = list(args)	# necessary bc of pop operation below
		
		# check that requested identifiers are present
		identifiers = self.tags.keys()
		for (i,key) in enumerate(args):
			if key not in identifiers:
				print 'WARNING: Ignoring ' + key + ' because it is not a recognized identifier.'
				args.pop(i)
				
		idxs = set([])
		for key in args:
			idxs.update(self.tags[key])
		idxs = list(idxs)
		
		return idxs
	
	def get_idxs_fullVDJ(self):
		"""Return indexes of chains that have a V, D, and J alignment."""
		# indexing of reference lists starts at 1 bc the first elt is ''
		Vs = set(self.get_idxs_OR(refseq.IGHV_list[1:]))
		Ds = set(self.get_idxs_OR(refseq.IGHD_list[1:]))
		Js = set(self.get_idxs_OR(refseq.IGHJ_list[1:]))
		idxs = list(Vs & Ds & Js)
		idxs.sort()
		return idxs
	
	def get_idxs_fullVJ(self):
		"""Return indexes of chains that have a V and J alignment."""
		Vs = set(self.get_idxs_OR(refseq.IGHV_list[1:]))
		Js = set(self.get_idxs_OR(refseq.IGHJ_list[1:]))
		idxs = list(Vs & Js)
		idxs.sort()
		return idxs
	
	def get_idxs_fullVJCDR3(self):
		"""Return indexes of chains that have a V and J alignment and a non-empty junction of length = 0 mod 3."""
		VJs   = set(self.get_idxs_fullVJ())
		CDR3s = set([i for (i,chain) in enumerate(self.chains) if chain.cdr3 % 3 == 0 and chain.cdr3 > 0])
		idxs = list(VJs & CDR3s)
		idxs.sort()
		return idxs		
	
	def get_chains_fullVDJ(self):
		"""Return subRepertoire that matches chains that have a V, D, and J alignment."""
		return self[ self.get_idxs_fullVDJ() ]
	
	def get_chains_fullVJ(self):
		"""Return subRepertoire that matches chains that have a V and J alignment."""
		return self[ self.get_idxs_fullVJ() ]
	
	def get_chains_fullVJCDR3(self):
		"""Return subRepertoire that matches chains that have a V and J alignment and a non-empty junction of length = 0 mod 3."""
		return self[ self.get_idxs_fullVJCDR3() ]
	
	# ==================================
	# = Appending/extending repertoire =
	# ==================================
	
	def __add__(self,other):
			return Repertoire(np.append(self.chains,other.chains), (self.metatags | other.metatags) )
	
	def __radd__(self,other):
		return Repertoire.__add__(other,self)
	
	def __iadd__(self,other):
		self.__init__( np.append(self.chains,other.chains), (self.metatags | other.metatags) )
		return self
	
	def __len__(self):
		return len(self.chains)
	
	def append(self,chain):
		print "WARNING: This function uses Repertoire.append(), which is slow and memory intensive for large repertoires."
		i = len(self.chains)
		self.chains = np.append(self.chains,chain)	# slow and inefficient: np.append is NOT in-place; it reallocates everytime
		self.processTags(i,chain)
		
		# not necessary to uniqueify the tags as each
		# time I add another one it's a higher unrepresented
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
		if isinstance(tagset,types.StringTypes): tagset = [tagset]
		for (i,chain) in enumerate(self.chains):
			chain.tags.update(tagset)
			for tag in tagset:
				try: self.tags[tag] += [i]
				except KeyError,e: self.tags[tag] = [i]
		self.uniqueifyTags()
		return
	
	def del_tags(self,tagset):
		if isinstance(tagset,types.StringTypes): tagset = [tagset]
		for tag in tagset:
			if tag not in self.tags.keys():
				print 'WARNING: ' + tag + 'is not a recognized tag'
				continue
			idxs = self.tags[tag]
			print idxs
			for chain in self.chains[idxs]:
				chain.tags.remove(tag)
			del self.tags[tag]
		return
	
	def add_metatags(self,metatagset):
		if isinstance(metatagset,types.StringTypes): metatagset = [metatagset]
		self.metatags.update(metatagset)
		return
	
	def del_metatags(self,metatagset):
		if isinstance(metatagset,types.StringTypes): metatagset = [metatagset]
		for tag in metatagset:
			if tag not in self.metatags:
				print 'WARNING: ' + tag + 'is not a recognized tag'
				continue
			self.metatags.remove(tag)
		return
		
	def processTags(self,i,chain):
		"""Register ImmuneChain tags in Repertoire object.
		
		Given ImmuneChain and position in array, get tags and add
		them to the tags dict of self
		
		"""
		
		for tag in chain.tags:
			try: self.tags[tag] += [i]
			except KeyError,e: self.tags[tag] = [i]
		
		try: self.tags[chain.v] += [i]
		except KeyError,e: self.tags[chain.v] = [i]
		
		try: self.tags[chain.d] += [i]
		except KeyError,e: self.tags[chain.d] = [i]
		
		try: self.tags[chain.j] += [i]
		except KeyError,e: self.tags[chain.j] = [i]
		
		try: self.tags[chain.ighc] += [i]
		except KeyError,e: self.tags[chain.ighc] = [i]
		
		try: self.tags[chain.descr] += [i]
		except KeyError,e: self.tags[chain.descr] = [i]
		
		return
	
	def uniqueifyTags(self):
		"""Process tag dict so there are no repeat indices."""
		for tag in self.tags.keys():
			self.tags[tag] = list(set(self.tags[tag]))
		return
	
	def reprocessTags(self):
		"""Process and uniqueify all tags of a Repertoire object from scratch."""
		self.tags = {}
		for ident in refseq.ALL_IDs[1:]:
			self.tags[ident] = []
		for (i,chain) in enumerate(self.chains):
			self.processTags(i,chain)
		self.uniqueifyTags()
		return
	
	def __str__(self):
		return self.__repr__()
	
	def __repr__(self):
		return self.getXML()
	
	def getXML(self,verbose=True):
		print "WARNING: the getXML function is slow and memory intensive for large repertoires."
		xmlstring = ''
		xmlstring += '<Repertoire>\n\n'
		xmlstring += '<Meta>\n'
		for tag in self.metatags:
			xmlstring += '\t<metatag>' + tag + '</metatag>\n'
		xmlstring += '</Meta>\n\n'
		xmlstring += '<Data>\n'
		for (i,chain) in enumerate(self.chains):
			if verbose and i%5000==0:
				print "Writing: " + str(i)
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
		self.saveData = False
		
		self.chains = []
		self.metatags = []
		
		self.chain_elements = [
							'descr',
							'seq',
							'v',
							'd',
							'j',
							'ighc',
							'cdr3',
							'junction',
							'func',
							'tag'
							  ]
	
	def startElement(self,name,attrs):
		
		#DEBUG
		#print "Start:  " + name
		#print "Buffer: " + '|'+self.buffer+'|'
		
		if name == 'ImmuneChain':
			self.currentChain = ImmuneChain()
		elif name == 'metatag':
			#DEBUG
			#print "started metatag"
			#print "\tbuffer: " + '|' + self.buffer + '|'
			self.saveData = True
		elif name in self.chain_elements:
			#DEBUG print 'caught elt: ' + name
			self.saveData = True
	
	def characters(self,content):
		if self.saveData:
			self.buffer += content
	
	def endElement(self,name):
		if name == 'metatag':
			#DEBUG
			#print "end metatag"
			#print "\tbuffer: " + '|' + self.buffer + '|'
			if self.mode == 'Repertoire':
				self.metatags.append(self.buffer)
		elif name == 'ImmuneChain':
			self.chains.append(self.currentChain)
			self.numChains += 1
		elif name in self.chain_elements:
			if name == 'cdr3':
				self.currentChain.cdr3 = eval(self.buffer)
			elif name == 'tag':
				self.currentChain.add_tags(self.buffer)
			else:
				self.currentChain.__setattr__(name,self.buffer)
		else:
			self.buffer = ''
		self.buffer = ''
		self.saveData = False
	
	def endDocument(self):
		print "READING VDJ XML"
		if self.mode == 'Repertoire':
			self.data.__init__(self.chains,self.metatags)
			print "mode: Repertoire"
			print "metatags:"
			for tag in self.data.metatags:
				print '\t' + tag
		elif self.mode == 'ImmuneChain':
			self.data.extend(self.chains)
			print "mode: ImmuneChain"
		print "Number of ImmuneChain objects read: " + str(self.numChains) + '\n'

def readVDJ(inputfile,mode='Repertoire'):
	"""Load a data from a VDJXML file as a Repertoire or list of ImmuneChains
	
	NOTE: this is MUCH slower than fastreadVDJ, which doesn't utilize
	the XML library itself.  This should only be used for more finicky files.
	
	"""
	
	if isinstance(inputfile,types.StringTypes):
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
	
	if isinstance(inputfile,types.StringTypes):
		ip.close()
	
	return data

def fastreadVDJ(inputfile,mode='Repertoire',verbose=True):
	"""Load a data from a VDJXML file as a Repertoire or list of ImmuneChains
	
	NOTE: this fn does NOT utilize the XML libraries; it implements a manual parser
	that takes input line by line.
	
	THIS ASSUMES THAT EVERY XML ELEMENT TAKES ONE AND ONLY ONE LINE
	
	"""
	if isinstance(inputfile,types.StringTypes):
		ip = open(inputfile,'r')
	elif isinstance(inputfile,file):
		ip = inputfile
	
	if mode == 'Repertoire':
		metatags = []
		data = []
	elif mode == 'ImmuneChain':
		data = []
	
	numChains = 0
	
	possible_elements = [
				'descr',
				'seq',
				'v',
				'd',
				'j',
				'ighc',
				'cdr3',
				'junction',
				'func',
				'tag'
				]
	
	for line in ip:
		line = line.strip()
		endelementpos = line.find('>') + 1
		xmlelement = line[0:endelementpos]
		element = xmlelement[1:-1]
		
		if xmlelement == '<metatag>':
			if mode == 'Repertoire':
				metatags.append(line[endelementpos:-1*(endelementpos+1)])
		elif xmlelement == '<ImmuneChain>':
			chain = ImmuneChain()
		elif xmlelement == '</ImmuneChain>':
			data.append(chain)
			numChains += 1
		elif element in possible_elements:
			if element == 'cdr3':
				chain.cdr3 = eval(line[endelementpos:-1*(endelementpos+1)])
			elif element == 'tag':
				chain.add_tags(line[endelementpos:-1*(endelementpos+1)])
			else:
				chain.__setattr__(element,line[endelementpos:-1*(endelementpos+1)])
	
	if isinstance(inputfile,types.StringTypes):
		ip.close()
	
	if mode == 'Repertoire':
		rep = Repertoire(data,metatags)
	elif mode == 'ImmuneChain':
		rep = data
	
	if verbose == True:
		print "READING VDJ XML"
		if mode == 'Repertoire':
			print "mode: Repertoire"
			print "metatags:"
			for tag in rep.metatags:
				print '\t' + tag
		elif mode == 'ImmuneChain':
			print "mode: ImmuneChain"
		print "Number of ImmuneChain objects read: " + str(numChains) + '\n'
	
	return rep

def writeVDJ(data, fileobj, verbose=True):
	"""Write Repertoire or list of ImmuneChains to a file in VDJXML format"""
	if isinstance(fileobj,types.StringTypes):
		handle = open(fileobj,'w')
	elif isinstance(fileobj,file):
		handle = fileobj
	
	if verbose == True:
		print "WRITING VDJ XML"
		if isinstance(data,Repertoire):
			print "mode: Repertoire"
			print "metatags:"
			for tag in data.metatags:
				print '\t' + tag
		elif isinstance(data,list) and isinstance(data[0],ImmuneChain):
			print "mode: ImmuneChain"
	
	print >>handle, '<?xml version="1.0"?>'
	if isinstance(data,Repertoire):
		# header
		xmlstring = ''
		xmlstring += '<Repertoire>\n\n'
		xmlstring += '<Meta>\n'
		for tag in data.metatags:
			xmlstring += '\t<metatag>' + tag + '</metatag>\n'
		xmlstring += '</Meta>\n\n'
		xmlstring += '<Data>\n'
		print >>handle, xmlstring
		
		# data
		for (i,chain) in enumerate(data.chains):
			#if verbose and i%5000==0:
			#	print "Writing: " + str(i)
			print >>handle, str(chain) + '\n'
		xmlstring = '</Data>\n\n'
		xmlstring += '</Repertoire>\n'
		print >>handle, xmlstring
	elif isinstance(data,list) and isinstance(data[0],ImmuneChain):
		for (i,chain) in enumerate(data):
			#if verbose and i%5000==0:
			#	print "Writing: " + str(i)
			print >>handle, chain
	
	if verbose == True:
		print "Number of ImmuneChain objects written: " + str(len(data)) + '\n'
	
	if isinstance(fileobj,types.StringTypes):
		handle.close()


############################################################################
############################################################################
###
###     CLEAN UP TO HERE
###
############################################################################
############################################################################




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

def countsClusters(clusters,reference_clusters):
	"""Takes a dictionary of cluster names mapped to a sequence of indexes in a repertoire.
	
	Returns an np array of the same length as reference_clusters with the counts of each
	cluster in reference_clusters.
	
	The need for reference_clusters is due to the fact that splitting a given repertoire
	may result in some parts not observing any clusters, so there needs to be a common way
	to compare two cluster sets
	
	"""
	
	counts = np.zeros(len(reference_clusters))
	for (i,name) in enumerate(reference_clusters):
		counts[i] = len(clusters.get(name,[]))
	return counts





class vdj_aligner(object):
	
	def __init__(self, verbose=False):
		
		t0 = time.time()
		
		self.numCrudeVCandidates = 5
		self.numCrudeDCandidates = 10
		self.numCrudeJCandidates = 2
		
		# Define seed patterns
		patternA='111011001011010111'
		patternB='1111000100010011010111'
		patternC='111111111111'
		patternD='110100001100010101111'
		patternE='1110111010001111'
		self.seedpatterns = [patternA,patternB,patternC,patternD,patternE]
		self.miniseedpatterns = ['111011','110111']
		
		# Generate hashes from reference data
		self.Vseqlistannot,self.Vseqlistkeys = vdj_aligner.seqlist2kmerannot( refseq.IGHV_seqs, self.seedpatterns )
		self.Dseqlistannotmini,self.Dseqlistkeysmini = vdj_aligner.seqlist2kmerannot( refseq.IGHD_seqs, self.miniseedpatterns )
		self.Jseqlistannot,self.Jseqlistkeys = vdj_aligner.seqlist2kmerannot( refseq.IGHJ_seqs, self.seedpatterns )
		
		t1 = time.time()
		
		if verbose: print "Database init:", t1-t0
		
		return
	
	def align_rep(self,rep,verbose=False):
		return self.align_chains(rep,verbose)
	
	def align_chains(self,chains,verbose=False):
		"""Align list of ImmuneChain objects to VDJ reference."""
		
		tstart = time.time()
		chaintimes = []
		vdj_junc = []
		
		for chain in chains:
			curraln = self.align_chain(chain)
			vdj_junc.append(curraln)
			if verbose: chaintimes.append(time.time())
		
		tend = time.time()
		
		if verbose:
			chaindurs = [(chaintimes[i+1]-chaintimes[i]) for i in xrange(len(chaintimes)-1)]
			print "Number of chains:", len(chains)
			print "Total time for aln:", tend-tstart
			print "Average time per chain:", np.mean(chaindurs)
		
		return vdj_junc
	
	def align_chain(self,chain,verbose=False):
		
		query = chain.seq
		
		t0 = time.time()
		
		# compute hashes from query seq
		queryannot,querykeys = vdj_aligner.seq2kmerannot(query,self.seedpatterns)
		
		t1 = time.time()
		
		Vscores_hash = {}
		
		# for each reference V segment and each pattern, how many shared k-mers are there?
		for Vseg in refseq.IGHV_seqs.keys():
			score = 0
			for pattern in self.seedpatterns:
				score += len( self.Vseqlistkeys[Vseg][pattern] & querykeys[pattern] )
			Vscores_hash[Vseg] = score
		
		# get numCrudeVCandidates highest scores in Vscores and store their names in descending order
		goodVseglist = [ seg[0] for seg in vdj_aligner.dict2sorteddecreasingitemlist(Vscores_hash,'value')[0:self.numCrudeVCandidates] ]
		
		t2 = time.time()
		
		bestVseg = ''
		bestVscore = 0
		bestVscoremat = []
		bestVtracemat = []
		
		# perform Needleman-Wunsch on top V seg candidates and remember which had the highest score
		for goodVseg in goodVseglist:
			# C implementation:
			# carve out memory
			# note that we are using zero initial conditions, so matrices are initialize too
			# notation is like Durbin p.29
			seq1 = refseq.IGHV_seqs[goodVseg]
			seq2 = query
			scores  = np.zeros( [len(seq1)+1, len(seq2)+1] )
			Ix = np.zeros( [len(seq1)+1, len(seq2)+1] )
			Iy = np.zeros( [len(seq1)+1, len(seq2)+1] )
			trace = np.zeros( [len(seq1)+1, len(seq2)+1], dtype=np.int)
			alignmentcore.alignNW( scores, Ix, Iy, trace, seq1, seq2 )
			
			# pure python:
			# scores,trace = vdj_aligner.alignNW( refseq.IGHV_seqs[goodVseg], query )
			
			currscore = vdj_aligner.scoreVJalign(scores)
			if currscore > bestVscore:
				bestVscore = currscore
				bestVseg = goodVseg
				bestVscoremat = scores
				bestVtracemat = trace
		
		chain.v = bestVseg
		
		t3 = time.time()
		
		# reconstruct the alignment and chop off V region through beginning of CDR3 (IMGT)
		if bestVseg != '':
			Valnref,Valnrefcoords,Valnquery,Valnquerycoords = vdj_aligner.construct_alignment( refseq.IGHV_seqs[bestVseg], query, bestVscoremat, bestVtracemat )
			query = vdj_aligner.pruneVregion( Valnref, Valnrefcoords, Valnquery, Valnquerycoords, bestVseg, query )
		
		t4 = time.time()
		
		# compute hashes from (pruned) query seq (junction + J)
		queryannot,querykeys = vdj_aligner.seq2kmerannot(query,self.seedpatterns)
		
		t5 = time.time()
		
		Jscores_hash = {}
		
		# for each reference J segment and each pattern, how many shared k-mers are there?
		for Jseg in refseq.IGHJ_seqs.keys():
			score = 0
			for pattern in self.seedpatterns:
				score += len( self.Jseqlistkeys[Jseg][pattern] & querykeys[pattern] )
			Jscores_hash[Jseg] = score
		
		# get numCrudeJCandidates highest scores in Jscores and store their names in descending order
		goodJseglist = [ seg[0] for seg in vdj_aligner.dict2sorteddecreasingitemlist(Jscores_hash,'value')[0:self.numCrudeJCandidates] ]
		
		t6 = time.time()
		
		bestJseg = ''
		bestJscore = 0
		bestJscoremat = []
		bestJtracemat = []
		
		# perform Needleman-Wunsch on top J seg candidates and remember which had the highest score
		for goodJseg in goodJseglist:
			# C implementation:
			# carve out memory
			# note that we are using zero initial conditions, so matrices are initialize too
			# notation is like Durbin p.29
			seq1 = refseq.IGHJ_seqs[goodJseg]
			seq2 = query
			scores  = np.zeros( [len(seq1)+1, len(seq2)+1] )
			Ix = np.zeros( [len(seq1)+1, len(seq2)+1] )
			Iy = np.zeros( [len(seq1)+1, len(seq2)+1] )
			trace = np.zeros( [len(seq1)+1, len(seq2)+1], dtype=np.int)
			alignmentcore.alignNW( scores, Ix, Iy, trace, seq1, seq2 )
			
			# pure python:
			#scores,trace = vdj_aligner.alignNW( refseq.IGHJ_seqs[goodJseg], query )
			
			currscore = vdj_aligner.scoreVJalign(scores)
			if currscore > bestJscore:
				bestJscore = currscore
				bestJseg = goodJseg
				bestJscoremat = scores
				bestJtracemat = trace
		
		chain.j = bestJseg
		
		t7 = time.time()
		
		# reconstruct the alignment and chop off J region after the TRP (IMGT)
		if bestJseg != '':
			Jalnref,Jalnrefcoords,Jalnquery,Jalnquerycoords = vdj_aligner.construct_alignment( refseq.IGHJ_seqs[bestJseg], query, bestJscoremat, bestJtracemat )
			query = vdj_aligner.pruneJregion( Jalnref, Jalnrefcoords, Jalnquery, Jalnquerycoords, bestJseg, query )
		
		t8 = time.time()
		
		# only attempt D alignment if both V and J were successful
		if bestVseg != '' and bestJseg != '':
			chain.junction = query
			chain.cdr3 = len(query)
			
			# compute hashes from junction sequence using mini seed patterns
			queryannot,querykeys = vdj_aligner.seq2kmerannot(query,self.miniseedpatterns)
		
			t9 = time.time()
		
			Dscores_hash = {}
		
			# for each reference D segment and each pattern, how many shared k-mers are there?
			for Dseg in refseq.IGHD_seqs.keys():
				score = 0
				for pattern in self.miniseedpatterns:
					score += len( self.Dseqlistkeysmini[Dseg][pattern] & querykeys[pattern] )
				Dscores_hash[Dseg] = score
		
			# get numCrudeDCandidates highest scores in Dscores and store their names in descending order
			goodDseglist = [ seg[0] for seg in vdj_aligner.dict2sorteddecreasingitemlist(Dscores_hash,'value')[0:self.numCrudeDCandidates] ]
		
			t10 = time.time()
		
			bestDseg = ''
			bestDscore = 0
			bestDscoremat = []
			bestDtracemat = []
		
			# perform Smith-Waterman on top D seg candidates and remember which had the highest score
			for goodDseg in goodDseglist:
				# C implementation:
				# carve out memory
				# note that we are using zero initial conditions, so matrices are initialize too
				# notation is like Durbin p.29
				seq1 = refseq.IGHD_seqs[goodDseg]
				seq2 = query
				scores  = np.zeros( [len(seq1)+1, len(seq2)+1] )
				trace = np.zeros( [len(seq1)+1, len(seq2)+1], dtype=np.int)
				alignmentcore.alignSW( scores, trace, seq1, seq2 )
				
				# pure python:
				#scores,trace = vdj_aligner.alignSW( refseq.IGHD_seqs[goodDseg], query )
				
				currscore = vdj_aligner.scoreDalign(scores)
				if currscore > bestDscore:
					bestDscore = currscore
					bestDseg = goodDseg
					bestDscoremat = scores
					bestDtracemat = trace
		
			t11 = time.time()
		
			chain.d = bestDseg
			
		else:
			bestDseg = ''
			t9 = t8
			t10 = t8
			t11 = t8
		
		# TODO: Rather than just assigning the best scorer to the V, D, or J segment
		#		I need a way to determine the alignment's significance, and only make
		#		make the assignment if there is a good match.
		
		if verbose:
			print t1-t0, "Full query hashes"
			print t2-t1, "Comparing hashes to V database"
			print t3-t2, "V seg NW alignment"
			print t4-t3, "Construct alignment and prune V region off"
			print t5-t4, "Pruned query hashes"
			print t6-t5, "Comparing hashes to J database"
			print t7-t6, "J seg NW alignment"
			print t8-t7, "Construct alignment and prune J region off"
			print t9-t8, "Pruned query hashes (junction only)"
			print t10-t9, "Comparing hashes to D database"
			print t11-t10, "D seg SW alignment"
			print t11-t0, "Total time"
		
		return chain.v, chain.d, chain.j, chain.junction
	
	@staticmethod
	def seq2kmerannot(seq1,patterns):
		"""Given sequence and patterns, for each pattern, compute all corresponding k-mers from sequence.
		
		The result is seqannot[pattern][key]=[pos1,pos2,...,posN] in seq
		
		"""
		seq = seqtools.seqString(seq1)
		seqannot = {}
		for pattern in patterns:
			seqannot[pattern] = {}
		
		for i in xrange(len(seq)):
			for pattern in patterns:
				word = seq[i:i+len(pattern)]
				if len(word) == len(pattern):
					key = ''.join( [p[1] for p in zip(pattern,word) if p[0]=='1'] )
				try: seqannot[pattern][key] += [i]
				except KeyError,e: seqannot[pattern][key] = [i]
		
		seqkeys = {}
		for pattern in patterns:
			seqkeys[pattern] = set( seqannot[pattern].keys() )
		
		return seqannot,seqkeys
	
	@staticmethod
	def seqlist2kmerannot(seqdict,patterns):
		seqlistannot = {}
		seqlistkeys  = {}
		for seq in seqdict.iteritems():
			seqannot,seqkeys = vdj_aligner.seq2kmerannot(seq[1],patterns)
			seqlistannot[seq[0]] = {}
			seqlistkeys[seq[0]]  = {}
			for pattern in patterns:
				seqlistannot[seq[0]][pattern] = seqannot[pattern]
				seqlistkeys[seq[0]][pattern]  = seqkeys[pattern]
		return seqlistannot,seqlistkeys
	
	@staticmethod
	def pruneVregion( alnref, alnrefcoords, alnquery, alnquerycoords, refID, queryseq ):
		"""Prune V region out of query sequence based on alignment.
		
		Given ref and query alignments of V region, refID, and the original
		query sequence, return a sequence with the V region cut out, leaving
		the 2nd-CYS.  Also needs query alignment coords.
		
		"""
		
		#DEBUG
		# # check that alnref actually has the whole reference segment
		# 		# otherwise, I would need to pass something like alnrefcoords
		# 		if alnref.replace('-','') != refseq.IGHV_seqs[refID]:
		# 			raise Exception, "Aligned reference segment is not equal to vdj.refseq reference segment."
		
		FR3end = refseq.IGHV_offset[refID] - alnrefcoords[0]		# first candidate position	
		#FR3end = refseq.IGHV_offset[refID]		# first candidate position	
		refgaps = alnref[:FR3end].count('-')	# count gaps up to putative CYS pos
		seengaps = 0
		while refgaps != 0:		# iteratively find all gaps up to the CYS
			seengaps += refgaps
			FR3end   += refgaps		# adjust if for gaps in ref alignment
			refgaps   = alnref[:FR3end].count('-') - seengaps	# any add'l gaps?
		
		querygaps = alnquery[:FR3end].count('-')
		
		# v_end_idx = idx of start of aln of query + distance into aln - # of gaps
		v_end_idx = alnquerycoords[0] + FR3end - querygaps
		
		return queryseq[v_end_idx:]
	
	@staticmethod
	def pruneJregion( alnref, alnrefcoords, alnquery, alnquerycoords, refID, queryseq ):
		"""Prune J region out of query sequence based on alignment.
		
		Given ref and query alignments of J region, refID, and the original
		query sequence, return a sequence with the J region cut out, leaving
		the J-TRP.  Also needs query alignment coords.
		
		"""
		#DEBUG
		# # check that alnref actually has the whole reference segment
		# 		# otherwise, I would need to pass something like alnrefcoords
		# 		if alnref.replace('-','') != refseq.IGHJ_seqs[refID]:
		# 			raise Exception, "Aligned reference segment is not equal to vdj.refseq reference segment."
		
		FR4start = refseq.IGHJ_offset[refID] - alnrefcoords[0]	# first candidate position of J-TRP start	
		refgaps = alnref[:FR4start].count('-')	# count gaps up to putative TRP pos
		seengaps = 0
		while refgaps != 0:		# iteratively find all gaps up to the TRP
			seengaps += refgaps
			FR4start += refgaps		# adjust for gaps in ref alignment
			refgaps   = alnref[:FR4start].count('-') - seengaps	# any add'l gaps?
		
		querygaps = alnquery[:FR4start].count('-')
		
		# v_end_idx = idx of start of aln of query + distance into aln - # of gaps + 3 nt for J-TRP
		j_start_idx = alnquerycoords[0] + FR4start - querygaps + 3
		
		return queryseq[:j_start_idx]
	
	@staticmethod
	def construct_alignment(seq1,seq2,scoremat,tracemat):
		"""Construct alignment of ref segment to query from score and trace matrices."""
		nrows,ncols = scoremat.shape
		
		# do some error checking
		if len(seq1)+1 != nrows or len(seq2)+1 != ncols:
			raise Exception, "nrows and ncols must be equal to len(seq1)+1 and len(seq2)+1"
		
		#DEBUG
		# if not nrows <= ncols:
		# 			raise Exception, "score matrix must have nrows < ncols"
		# 		if not len(seq1) <= len(seq2):
		# 			raise Exception, "len of seq1 must be smaller than seq2"
		
		# translate integer traces to coords
		deltas = {
			0 : (1,1),
			1 : (1,0),
			2 : (0,1),
			3 : (0,0)
		}
		
		# compute col where alignment should start
		if nrows <= ncols:
			col = np.argmax( scoremat[nrows-1,:] )
			row = nrows-1
		else:
			col = ncols-1
			row = np.argmax( scoremat[:,ncols-1] )
		#DEBUG
		# col = np.argmax( scoremat[nrows-1,:] )
		# row = nrows-1
		
		# if row is coord in matrix, row-1 is coord in seq
		
		aln1 = seq1[row-1]
		aln2 = seq2[col-1]
		
		aln1end = row
		aln2end = col
		
		#DEBUG
		#while row-1 > 0:
		while (row-1 > 0) and (col-1 > 0):
			# compute direction of moves
			rowchange,colchange = deltas[ tracemat[row,col] ]
			
			# WORKS WITH PURE PYTHON alignNW trace return
			#rowchange = row-tracemat[row,col][0]
			#colchange = col-tracemat[row,col][1]
			
			# emit appropriate symbols
			if rowchange == 1:
				row -= 1
				aln1 = seq1[row-1] + aln1
			elif rowchange == 0:
				aln1 = '-' + aln1
			else:
				raise Exception, "Trace matrix contained jump of greater than one row/col."
			
			if colchange == 1:
				col -= 1
				aln2 = seq2[col-1] + aln2
			elif colchange == 0:
				aln2 = '-' + aln2
			else:
				raise Exception, "Trace matrix contained jump of greater than one row/col."
		
		aln1start = row-1
		aln2start = col-1
		
		return aln1, (aln1start,aln1end), aln2, (aln2start,aln2end)
	
	@staticmethod
	def scoreVJalign(scorematrix):
		"""Computes score of V alignment given Needleman-Wunsch score matrix
		
		ASSUMES num rows < num cols, i.e., refseq V seg is on vertical axis
		
		"""
		nrows,ncols = scorematrix.shape
		
		
		if nrows <= ncols:
			return np.max( scorematrix[nrows-1,:] )
		else:
			return np.max( scorematrix[:,ncols-1] )
		
		#DEBUG
		#OLD WAY
		# nrows,ncols = scorematrix.shape
		# 		
		# 		if not nrows < ncols:
		# 			raise Exception, "score matrix must have nrows < ncols"
		# 		
		# 		return np.max( scorematrix[nrows-1,:] )
		
	@staticmethod
	def scoreDalign(scorematrix):
		"""Computes score of D alignment given Smith-Waterman score matrix
		
		"""
		return np.max( scorematrix )
	
	@staticmethod
	def alignNW(sequence1,sequence2):
		'''Align two seqs using modified Needleman-Wunsch algorithm as in Durbin.
	
		Returns only a single optimal alignment in the backtrace.
		
		Initialization is with zeros rather than with gap penalties.  This allows
		the smaller reference segment to freely align anywhere along the query
		sequence.
	
		'''
		
		seq1 = seqtools.seqString(sequence1)
		seq2 = seqtools.seqString(sequence2)
		
		# define parameters
		match     =  0.5
		mismatch  = -0.75
		gapopen   = -2.0
		gapextend = -1.5
		
		def match_fn(a,b):
			if a == b:
				return match
			else:
				return mismatch
		
		# carve out memory
		# note that we are using zero initial conditions, so matrices are initialize too
		# notation is like Durbin p.29
		M  = np.zeros( [len(seq1)+1, len(seq2)+1] )
		Ix = np.zeros( [len(seq1)+1, len(seq2)+1] )
		Iy = np.zeros( [len(seq1)+1, len(seq2)+1] )
		BT = np.zeros( [len(seq1)+1, len(seq2)+1], dtype=np.object)
		
		#alignmentcore.alignNW(M,Ix,Iy,BT,seq1,seq2)
		
		for i in xrange( 1, len(seq1)+1 ):
			for j in xrange( 1, len(seq2)+1 ):
				s = match_fn(seq1[i-1],seq2[j-1])	# note that seqs are zero-indexe
				extensions = [ M[i-1,j-1] + s, Ix[i-1,j-1] + s, Iy[i-1,j-1] + s ]	# Durbin
				#extensions = [ M[i-1,j-1] + s, Ix[i-1,j-1], Iy[i-1,j-1] ]	# Gotoh
				pointers = [ (i-1,j-1), (i-1,j), (i,j-1) ]
				best = np.argmax( extensions )
			
				M[i,j]  = extensions[best]
				Ix[i,j] = max( M[i-1,j] + gapopen, Ix[i-1,j] + gapextend )
				Iy[i,j] = max( M[i,j-1] + gapopen, Iy[i,j-1] + gapextend )
			
				BT[i,j] = pointers[ best ]
		
		return M, BT
	
	@staticmethod
	def alignSW(sequence1,sequence2):
		"""Align two seqs using classic Smith-Waterman algorithm as in Durbin.
		
		Uses linear gap costs
		
		Returns only a single optimal alignment in the backtrace.
		
		Notation is like the original paper.
		
		"""
		
		seq1 = seqtools.seqString(sequence1)
		seq2 = seqtools.seqString(sequence2)
		
		# define parameters
		match     =  0.5
		mismatch  = -0.75
		gapextend = -1.5
		
		def match_fn(a,b):
			if a == b:
				return match
			else:
				return mismatch
		
		# carve out memory
		# note that we are using zero initial conditions, so matrices are initialized too
		# notation is like Durbin p.22
		F = np.zeros( [len(seq1)+1, len(seq2)+1] )
		BT = np.zeros( [len(seq1)+1, len(seq2)+1], dtype=np.object)
		
		for i in xrange( 1, len(seq1)+1 ):
			for j in xrange( 1, len(seq2)+1 ):
				s = match_fn(seq1[i-1],seq2[j-1])	# note that seqs are zero-indexe
				extensions = [ 0, F[i-1,j-1] + s, F[i-1,j] + gapextend, F[i,j-1] + gapextend ]
				pointers = [ 0, (i-1,j-1), (i-1,j), (i,j-1) ]
				best = np.argmax( extensions )
			
				F[i,j]  = extensions[ best ]
				BT[i,j] = pointers[ best ]
		
		return F, BT
	
	@staticmethod
	def dict2sorteddecreasingitemlist(dictionary,keyorvalue='value'):
		pos = {'key':0, 'value':1}
		di = dictionary.items()
		di.sort(key=operator.itemgetter(pos[keyorvalue]))
		di.reverse()
		return di
















# 
# 
# 
# 
# def alignV(query,Vseg):
# 	'''Align query sequence using modified Smith-Waterman alignment to reference V segment.
# 	
# 	Uses Smith-Waterman like algorithm with affine gap penalties.
# 	
# 	The alignment is ASYMMETRIC.  Overhang penalties are only in one direction.  Because the V seg
# 	should align to the 5' end of query, overhang on the left is penalized while overhang on the
# 	right is not.
# 	
# 	'''
# 	matchscore     =  0.5
# 	mismatchscore  = -0.75
# 	gapopenscore   = -2.0
# 	gapextendscore = -1.5
# 	
# 	scores    = np.zeros( [len(query)+1, len(Vseg)+1] )
# 	backtrace = np.zeros( [len(query)+1, len(Vseg)+1] )
# 	
# 	for r in xrange( len(query)+1 ):
# 		pass
# 	
# 	











# class abacus_aligner(object):
# 	def __init__(self,verbose=False):
# 		
# 		t0 = time.time()
# 		
# 		# Define seed patterns
# 		patternA='111011001011010111'
# 		patternB='1111000100010011010111'
# 		patternC='111111111111'
# 		patternD='110100001100010101111'
# 		patternE='1110111010001111'
# 		self.seedpatterns = [patternA,patternB,patternC,patternD,patternE]
# 		self.miniseedpatterns = ['111011','110111']
# 		
# 		# Load reference germline library
# 		self.Vrefseqlist = seqtools.getFasta(V_ref_fasta_file)
# 		self.Drefseqlist = seqtools.getFasta(D_ref_fasta_file)
# 		self.Jrefseqlist = seqtools.getFasta(J_ref_fasta_file)
# 		
# 		self.Vrefseqdict = seqlist2seqdict(Vrefseqlist)
# 		self.Drefseqdict = seqlist2seqdict(Drefseqlist)
# 		self.Jrefseqdict = seqlist2seqdict(Jrefseqlist)
# 		
# 		# Generate hashes from reference data
# 		self.Vseqlistannot,self.Vseqlistkeys = seqlist2kmerannot( self.Vrefseqlist, self.seedpatterns )
# 		self.Dseqlistannot,self.Dseqlistkeys = seqlist2kmerannot( self.Drefseqlist, self.seedpatterns )
# 		self.Jseqlistannot,self.Jseqlistkeys = seqlist2kmerannot( self.Jrefseqlist, self.seedpatterns )
# 		
# 		self.Dseqlistannotmini,self.Dseqlistkeysmini = seqlist2kmerannot( self.Drefseqlist, self.miniseedpatterns )
# 		
# 		t1 = time.time()
# 		
# 		if verbose:
# 			print "Database init:", t1-t0
# 		
# 		return
# 	
# 	# =============
# 	# = Interface =
# 	# =============
# 	
# 	def align_seq(self,seq,verbose=False):
# 		
# 		query = seqtools.seqString(seq)
# 		
# 		t0 = time.time()
# 		
# 		queryannot,querykeys = seq2kmerannot(query,self.seedpatterns)
# 		
# 		t1 = time.time()
# 		
# 		Vscores = {}
# 		Dscores = {}
# 		Jscores = {}
# 		
# 		# for each reference segment and each pattern, how many shared k-mers are there?
# 		for Vseg in self.Vrefseqdict.keys():
# 			score = 0
# 			for pattern in self.seedpatterns:
# 				score += len( self.Vseqlistkeys[Vseg][pattern] & querykeys[pattern] )
# 			Vscores[Vseg] = scores
# 		
# 		for Jseg in self.Jrefseqdict.keys():
# 			score = 0
# 			for pattern in self.seedpatterns:
# 				score += len( self.Jseqlistkeys[Jseg][pattern] & querykeys[pattern] )
# 			Jscores[Jseg] = scores
# 		
# 		goodVscores = {}
# 		goodJscores = {}
# 		bestV = ''
# 		bestJ = ''
# 		
# 		# get 10 highest scores in Vscores and store their names in descending order
# 		goodVseglist = [ seg[0] for seg in dict2sorteditemlist(Vscores,'value').reverse()[0:10] ]
# 	
# 	def align_seqlist(self,seqs):
# 		pass
# 	
# 	def align_fasta(self,fastafile):
# 		pass
# 	
# 	# =====================
# 	# = Utility functions =
# 	# =====================
# 	
# 	@staticmethod
# 	def seqlist2seqdict(seqlist):
# 		seqdict = {}
# 		for seq in seqlist:
# 			seqdict[seq.description] = seqtools.seqString(seq)
# 		return seqdict
# 	
# 	@staticmethod
# 	def seq2kmerannot(seq1,patterns):
# 		"""Given sequence and patterns, for each pattern, compute all corresponding k-mers from sequence.
# 		
# 		The result is seqannot[pattern][key]=[pos1,pos2,...,posN] in seq
# 		
# 		"""
# 		seq = seqtools.seqString(seq1)
# 		seqannot = {}
# 		for pattern in patterns:
# 			seqannot[pattern] = {}
# 		
# 		for i in xrange(len(seq)):
# 			for pattern in patterns:
# 				word = seq[i:i+len(pattern)]
# 				if len(word) == len(pattern):
# 					key = ''.join( [p[1] for p in zip(pattern,word) if p[0]=='1'] )
# 				try: seqannot[pattern][key] += [i]
# 				except KeyError,e: seqannot[pattern][key] = [i]
# 		
# 		seqkeys = {}
# 		for pattern in patterns:
# 			seqkeys[pattern] = set( seqannot[pattern].keys() )
# 		
# 		return seqannot,seqkeys
# 	
# 	@staticmethod
# 	def seqlist2kmerannot(seqlist,patterns):
# 		seqlistannot = {}
# 		seqlistkeys  = {}
# 		for seq in seqlist:
# 			seqannot,seqkeys = seq2kmerannot(seqtools.seqString(seq),patterns)
# 			seqlistannot[seq.description] = {}
# 			seqlistkeys[seq.description]  = {}
# 			for pattern in patterns:
# 				seqlistannot[seq.description][pattern] = seqannot[pattern]
# 				seqlistkeys[seq.description][pattern]  = seqkeys[pattern]
# 		return seqlistannot,seqlistkeys
# 	
# 	@staticmethod
# 	def dict2sorteditemlist(dictionary,keyorvalue='value'):
# 		pos = {'key':0, 'value':1}
# 		return dictionary.items().sort(key=operator.itemgetter[pos[keyorvalue]])