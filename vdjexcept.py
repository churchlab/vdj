# vdjexcept.py

import exceptions

class AlignmentError(exceptions.Exception):
	def __init__(self,msg,data):
		self.args = (msg,data)
		self.msg = msg
		self.data = data


class NoData(AlignmentError):
	def __init__(self,msg,data):
		self.args = (msg,data)
		self.msg = msg
		self.data = data


class IncompleteData(AlignmentError):
	def __init__(self,msg,chain):
		self.args = (msg,chain)
		self.msg = msg
		self.data = chain
	


class NonReference(AlignmentError):
	def __init__(self,msg,chain):
		self.args = (msg,chain)
		self.msg = msg
		self.data = chain
	
