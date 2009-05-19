# __init__.py

# VDJ package: toolset for vdj manipulations

import numpy as np
import seqtools
import refseq
import clones	#CDR3 extraction
import alignment
import types
import math
import tempfile
import subprocess
import time
import os
import xml.sax
import xml.sax.handler
from xml.sax.saxutils import escape, unescape
from Bio.Seq import reverse_complement
from datetime import datetime
from editdist import distance as editdist
import scipy.cluster

NV = len(refseq.IGHV)
ND = len(refseq.IGHD)
NJ = len(refseq.IGHJ)

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
		if isinstance(tags,types.StringTypes): tags = [tags]
		self.tags = set(tags)	# tag for sample number/experiment etc
		self.v = v
		self.d = d
		self.j = j
		self.ighc = ighc
		self.cdr3 = cdr3
		self.junction = junction
		self.func = func
	
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
		
	def __init__( self, chains=[], metatags=[] ):
		'''
		must provide list of ImmuneChain objects
		they don't have to have full information
		'''
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
		for ident in refseq.IGHV + refseq.IGHD + refseq.IGHJ + refseq.IGHC:
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
		idxs = self.get_idxs_AND(args)
		return Repertoire(self.chains[idxs],self.metatags)
	
	def get_chains_OR(self,args):
		idxs = self.get_idxs_OR(args)		
		return Repertoire(self.chains[idxs],self.metatags)
		
	def get_idxs_AND(self,args):
		if isinstance(args,types.StringTypes): args = [args]
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
		if isinstance(args,types.StringTypes): args = [args]
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
		self.__init__( np.append(self.chains,other.chains), (self.metatags | other.metatags) )
		return self
	
	def __len__(self):
		return len(self.chains)
	
	def append(self,chain):
		print "WARNING: This function uses Repertoire.append(), which is slow and memory intensive for large repertoires."
		i = len(self.chains)
		self.chains = np.append(self.chains,chain)	# slow and inefficient: np.append is NOT in-place
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
		
		try: self.tags[chain.ighc] += [i]
		except KeyError,e: self.tags[chain.ighc] = [i]
		
		try: self.tags[chain.descr] += [i]
		except KeyError,e: self.tags[chain.descr] = [i]
		
		return
	
	def uniqueifyTags(self):
		for tag in self.tags.keys():
			self.tags[tag] = list(set(self.tags[tag]))
		return
	
	def reprocessTags(self):
		self.tags = {}
		for ident in refseq.IGHV + refseq.IGHD + refseq.IGHJ + refseq.IGHC:
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
	'''
	function to load a Repertoire and return it in a Repertoire object
	or a list of ImmuneChains
	'''
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
	'''
	function to load a Repertoire and return in a repertoire object
	it does not utilize the XML libraries, but manually reads the input
	file line by line.
	
	THIS ASSUMES THAT EVERY XML ELEMENT TAKES ONE AND ONLY LINE
	'''
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

#===============================================================================

# ======================
# = Pipeline functions =
# ======================

def initial_import(inputfilelist,outputname,metatags=[],tags=[]):
	seqs = []
	for filename in inputfilelist:
		seqs.extend(seqtools.getFasta(filename))
	
	chains = [ImmuneChain(descr=seq.description.split()[0],seq=seq.seq.data,tags=tags) for seq in seqs]	
	rep = Repertoire(chains)
	rep.add_metatags(metatags)
	rep.add_metatags("Initial_Import : " + timestamp())
	
	writeVDJ(rep,outputname)
	
	return rep

def size_select(rep,readlensizes):	
	if len(readlensizes) != 2:
		raise Exception, "Incorrect number of args for size_select operation."
	
	#DEBUG
	#print "size select: " + str(len(rep))
	
	minreadlen = readlensizes[0]
	maxreadlen = readlensizes[1]
	
	idxs = [i for (i,chain) in enumerate(rep) if len(chain) >= minreadlen and len(chain) <= maxreadlen]
	
	#DEBUG
	#print idxs
	#print str(type(idxs))
	
	rep_sizeselected = rep[idxs]
	rep_sizeselected.add_metatags("Size Selection: " +"min " + str(minreadlen)+ " max " + str(maxreadlen) + " : " + timestamp())
		
	return rep_sizeselected

def barcode_id(rep,barcodefile):
	# load the barcodes from file and process them
	ip = open(barcodefile,'r')
	barcodes = {}
	for line in ip:
		bc = line.split()
		barcodes[bc[0]] = bc[1]
	ip.close()
	barcodelen = len(barcodes.keys()[0])
	
	rep.add_metatags("Barcode_ID : " + timestamp())
	
	for chain in rep:
		curr_bcID = barcodes.get( chain.seq[:barcodelen], '' )
		if curr_bcID != '':
			chain.seq = chain.seq[barcodelen:]	# remove barcode
			chain.add_tags(curr_bcID)
	
	rep.reprocessTags()
	
	return rep

def isotype_id(rep,IGHCfile):
	# load the isotype ids from file and process them
	ip = open(IGHCfile,'r')
	isotypes = {}
	for line in ip:
		iso = line.split()
		# note that I take the revcomp of the primer sequence here
		isotypes[reverse_complement(iso[0])] = iso[1]
	ip.close()
	
	rep.add_metatags("Isotype_ID : " + timestamp())
	
	for chain in rep:
		curr_isoID = ''
		for iso in isotypes.iteritems():
			if iso[0] in chain.seq[-50:]:	# arbitrary cutoff from 3' end
				curr_isoID = iso[1]
		chain.ighc = curr_isoID
	
	rep.reprocessTags()
	
	return rep

def positive_strand(rep):
	rep.add_metatags("Positive_Strand : " + timestamp())
	
	# compute correct strand
	strands = alignment.chainlist2strandlist(rep)
	if len(strands) != len(rep):
		raise Exception, "Error in positive strand id."
	
	# invert strand
	for (strand,chain) in zip(strands,rep):
		chain.add_tags('positive')
		if strand == -1:
			chain.add_tags('revcomp')
			chain.seq = reverse_complement(chain.seq)
	
	rep.reprocessTags()
	
	return rep

def align_rep(rep):
	rep.add_metatags("VDJ_Alignment : " + timestamp())
	
	# perform alignment
	alignment.abacus.alignchainlist(rep)
	extractCDR3(rep)
	
	rep.reprocessTags()
	
	return rep

#===============================================================================

# ===================
# = LSF Dispatching =
# ===================

def split_into_parts(rep,outputname,packetsize):
	parts = []
	nfiles = int(math.ceil(float(len(rep)) / packetsize))
	for i in xrange(nfiles):
		# split into multiple parts; append .partNNN to filename
		currfilename = outputname + '.part%04i' % i
		parts.append(currfilename)
		writeVDJ(rep[i*packetsize:(i+1)*packetsize],currfilename)
	return parts

def load_parts(parts):
	rep = Repertoire()
	for part in parts:
		rep_part = fastreadVDJ(part,mode='Repertoire')
		rep += rep_part
	return rep

def generate_script(operation):
	scriptname = tempfile.mktemp('.py','vdj_operation','./')
	op = open(scriptname,'w')
	print >>op, "import vdj"
	print >>op, "import sys"
	print >>op, "rep = vdj.fastreadVDJ(sys.argv[1],mode='Repertoire')"
	print >>op, "rep=vdj."+operation+"(rep)"
	print >>op, "vdj.writeVDJ(rep,sys.argv[1])"
	op.close()
	return scriptname

def submit_to_LSF(queue,LSFopfile,script,parts):
	processes = []	# list of PIDs that are dispatched through LSF
	for chunk in parts:
		proc = subprocess.Popen( ['bsub','-q'+queue,'-o'+LSFopfile,'python',script,chunk], stdout=subprocess.PIPE )
		proc.wait()
		processes.append(proc.stdout.read().split('<')[1].split('>')[0])
	return processes

def waitforLSFjobs(PIDs,interval=30):
	finished = False
	while not finished:
		time.sleep(interval)
		p = subprocess.Popen('bjobs',stdout=subprocess.PIPE)
		p.wait()
		status = p.stdout.read().split('\n')
		if status[0].split()[0] != 'JOBID':
			finished = False
			continue
		runningprocesses = [line.split()[0] for line in status if line.split() != [] and line.split()[0] != 'JOBID']
		finished = True
		for pid in PIDs:
			if pid in runningprocesses:
				finished = False
	return

#===============================================================================

# ==================
# = Misc Utilities =
# ==================

def timestamp():
	return datetime.now().isoformat().split('.')[0]

#===============================================================================

# =============================
# = Clustering/CDR3 functions =
# =============================

# Distance metrics

def chain_Levenshtein(x,y):
	return editdist(x.seq,y.seq)

def chain_junction_Levenshtein(x,y):
	return editdist(x.junction,y.junction)

def NGLD(x,y):
	'''
	Normalized Generalized Levenshtein Distance
	(generalization is trivial case here; alpha=1)
	based on IEEE Trans Pattern Analys Mach Intel 29(6):1091
	'''
	try:
		GLD = editdist(x,y)
		Nx	= float(len(x))
		Ny	= float(len(y))
		return (2.*GLD)/(Nx+Ny+GLD)
	except ZeroDivisionError, e:
		return 0.

# Clustering functions

def pdist(X,metric):
	m = len(X)
	dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)
	k = 0
	for i in xrange(0, m - 1):
		for j in xrange(i+1, m):
			dm[k] = metric(X[i], X[j])
			k += 1
	return dm

def clusterChains(chains,cutoff=4.5,tag_chains=False,tag=''):
	Y = pdist(chains,chain_junction_Levenshtein)
	Z = scipy.cluster.hierarchy.linkage(Y,method='average')
	T = scipy.cluster.hierarchy.fcluster(Z,cutoff,criterion='distance')
	if tag_chains == True:
		for (i,chain) in enumerate(chains):
			chain.add_tags('cluster|'+tag+'|'+str(T[i]))
	return T

def clusterRepertoire(rep,cutoff=4.5,tag_chains=False,tag=''):
	'''
	split repertoire into all V-J combos and perform clustering
	tag_chains will add a cluster tag to each chain; tag will be
	incorporated as well.
	If false, then it will return a list of lists, each one of
	which represents a cluster and has the chain descr in it
	'''
	repgood = rep.get_chains_fullVJCDR3()
	clusters = []
	for vseg in refseq.IGHV[1:]:
		for jseg in refseq.IGHJ[1:]:
			currtag = tag+'|'+vseg+'|'+jseg
			currchains = repgood.get_chains_AND([vseg,jseg]).chains
			T = clusterChains(repgood.get_chains_AND([vseg,jseg]).chains,cutoff,tag_chains,currtag)
			numclusters = len(set(T))
			currclusters = [ [] for i in np.arange(numclusters)]
			for (i,clust) in enumerate(T):
				currclusters[clust].append(currchains[i].descr)
			clusters.extend(currclusters)
	return clusters

# CDR3 Extraction

def extractCDR3(chains):
	'''
	takes list of ImmuneChain objects.	They must have a sequence, a V, and J region alignment
	'''
	
	if isinstance(chains,Repertoire): chain_list = chains.chains
	if isinstance(chains,list): chain_list = chains
	if isinstance(chains,ImmuneChain): chain_list = [chains]
	
	curdir = os.getcwd()
	if curdir[-1] != '/': curdir += '/'	   
	
	query_file_name = 'query'+seqtools.randalphanum(8)+'.fasta'
	while os.path.exists(query_file_name):
		query_file_name = 'query'+seqtools.randalphanum(8)+'.fasta'
	
	ref_FR3_file_name = 'ref_FR3_'+seqtools.randalphanum(8)+'.fasta'
	while os.path.exists(ref_FR3_file_name):
		ref_FR3_file_name = 'ref_FR3_'+seqtools.randalphanum(8)+'.fasta'
	
	ref_J_file_name	  = 'ref_J_'+seqtools.randalphanum(8)+'.fasta'
	while os.path.exists(ref_J_file_name):
		ref_J_file_name	  = 'ref_J_'+seqtools.randalphanum(8)+'.fasta'
	
	cdr3s = []
	for chain in chain_list:
		try:
			# init indexes
			v_end_idx = len(chain.seq)
			j_start_idx = 0
			
			# can't get the CDR3
			if chain.v == '' or chain.j == '':
				cdr3s.append('')
				continue
			
			# write sequence and references to fasta file
			query_file	 = open(curdir+query_file_name, 'w')
			ref_FR3_file = open(curdir+ref_FR3_file_name,'w')
			ref_J_file	 = open(curdir+ref_J_file_name, 'w')
			
			print >>query_file, '>query_V'
			print >>query_file, chain.seq
			
			print >>ref_FR3_file, '>ref_V_FR3'
			print >>ref_FR3_file, refseq.IGHV_FR3[chain.v]
			
			print >>ref_J_file, '>ref_J_FR4'
			print >>ref_J_file, refseq.IGHJ_FR4[chain.j][1]
			
			query_file.close()
			ref_FR3_file.close()
			ref_J_file.close()
			
			v_align = os.popen('bl2seq -i ' + query_file_name + ' -j ' + ref_FR3_file_name + ' -p blastn -g T -F F -e 1.0 -D 1 -r 3 -G 3 -E 6 -W 9')
			
			for line in v_align:
				if line[0]=='#':
					continue
				linedata = line.split()
				
				refend   = eval(linedata[-3])
				reflen   = len(refseq.IGHV_FR3[chain.v])
				ext_deficit = reflen - refend
				
				# TODO: sometimes bl2seq will not extend to the end of FR3
				# to resolve this, I need to do the same thing I do with the
				# J region, and move the 2nd-CYS to an internal site.  the prob I
				# forsee here is that there is not much add'l conserved seq
				# after the 2nd-CYS
				if ext_deficit > 3:
					raise Exception, 'V segment FR3 alignment for sequence ' + chain.descr + ' was not extended through the end'
				queryend = eval(linedata[-5])
				v_end_idx = eval(line.split()[-5]) + ext_deficit - 3		# FR3 includes the CYS
				break
			
			query_file	 = open(curdir+query_file_name, 'w')
			
			print >>query_file, '>query_J'
			print >>query_file, chain.seq[v_end_idx:]
			
			query_file.close()
			
			j_align = os.popen('bl2seq -i ' + query_file_name + ' -j ' + ref_J_file_name + ' -p blastn -g T -F F -e 1.0 -D 0 -r 3 -G 3 -E 6 -W 9')
			
			# as the conserved TRP is internal to the ref sequence, we have to deal with possible gapped alignments
			for line in j_align:
				linedata = line.split()
				if linedata == [] or linedata[0] != 'Query:':
					continue
				
				querystart = eval(linedata[1])
				querymatch = linedata[2]
				
				j_align.next()	# burn a line
				line = j_align.next()
				linedata = line.split()
				
				refstart   = eval(linedata[1])
				refmatch   = linedata[2]
				
				# count the number of gaps before the position of the J-TRP
				TRPstart = refseq.IGHJ_FR4[chain.j][0] - (refstart - 1)	# first candidate pos
				if TRPstart < 0 or TRPstart > len(refmatch):
					raise Exception, 'J segment FR4 alignment for sequence ' + chain.descr + ' did not include the J-TRP site'
				refgaps  = refmatch[:TRPstart].count('-')	# count gaps up to putative TRP pos
				seengaps = 0
				while refgaps != 0:	# iteratively find all gaps up to TRP
					seengaps += refgaps
					TRPstart += refgaps	# adjust it for any gaps in ref seq
					refgaps   = refmatch[:TRPstart].count('-') - seengaps	# any add'l gaps?
				
				querygaps = querymatch[:TRPstart].count('-')	# num gaps in corres query segment
				
				# j_start_idx = pos i left off for V + distance into what's left + distance to TRP - adj for gaps + 3 for TRP
				j_start_idx = v_end_idx + (querystart - 1) + TRPstart - querygaps + 3
				
				break	# in case there are less high-scoring matches
			
			cdr3s.append(chain.seq[v_end_idx:j_start_idx])
		except Exception, e:
			print e
			cdr3s.append('')
			
	
	if os.path.exists(curdir+query_file_name): os.remove(curdir+query_file_name)
	if os.path.exists(curdir+ref_FR3_file_name): os.remove(curdir+ref_FR3_file_name)
	if os.path.exists(curdir+ref_J_file_name): os.remove(curdir+ref_J_file_name)
	
	if isinstance(chains,Repertoire):
		if len(chains.chains) != len(cdr3s):
			raise Exception, 'wrong number of CDR3s extracted'
		for (chain,cdr3) in zip(chains.chains,cdr3s):
			chain.junction = cdr3
			chain.cdr3 = len(cdr3)
	if isinstance(chains,list):
		if len(chains) != len(cdr3s):
			raise Exception, 'wrong number of CDR3s extracted'
		for (chain,cdr3) in zip(chains,cdr3s):
			chain.junction = cdr3
			chain.cdr3 = len(cdr3)
	if isinstance(chains,ImmuneChain):
		if len(cdr3s) != 1:
			raise Exception, 'wrong number of CDR3s extracted'
		chains.junction = cdr3s[0]
		chains.cdr3 = len(cdr3s[0])
	
	return chains

#===============================================================================


#===============================================================================
