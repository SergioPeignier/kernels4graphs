class atom:
	#name='C'#carbon ...
	#deg=2#number of neigboors
	#morgan=3 morgan coeficient	
	""" creat an atom object
		:param self: atom
		:param name : atom label 
		#further parameters can be added to this object
		"""
	def __init__(self,name):
		self.name=name
		self.deg=0
		self.morgan=1
		self.new_morgan=0
		self.old_name=name
	""" surcharge the print function
		:param self: atom
		"""
	def __str__(self):
		return "atom: "+self.name+" | deg: "+str(self.deg)+" | morgan: "+str(self.morgan)	
	""" surcharge the equality function
		:param self: another atom
		"""	
	def __eq__(self,other):
		return self.name==other.name
	""" surcharge the inequality function
		:param self: another atom
		"""
	def __ne__(self,other):
		return self.name!=other.name
	""" update the atom's morgan index
		:param self: atom 
		"""
	def refresh_morgan(self):
		self.morgan = self.new_morgan
		self.new_morgan=0
	""" include the morgan index into the atom's name
		:param self: atom
		"""
	def include_morgan_in_name(self):
		self.name=self.old_name+str(self.morgan)
