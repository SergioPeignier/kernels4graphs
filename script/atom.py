class atom:
	def __init__(self,name, deg=0, morgan=1):
		""" Create an atom object. name='C' for Carbon atom, deg holds for the number of neighbors, morgan morgan coeficient
		:param self: atom object
		:param name: atom label 
		:type name: String
		"""
		self.name=name
		self.deg=deg
		self.morgan=morgan
		self.new_morgan=0
		self.old_name=name
	
	def __str__(self):
		""" Redefine print function
		:param self: atom
		"""
		return "atom: "+self.name+" | deg: "+str(self.deg)+" | morgan: "+str(self.morgan)	
	
	def __eq__(self,other):
		""" Redefine the equality function
		:param self: this atom object
		:param other: another atom object
		:type other: atom
		"""	
		return self.name==other.name
	
	def __ne__(self,other):
		""" 
		Redefine the inequality function
		:param self: this atom object
		:param other: another atom object
		:type other: atom
		"""
		return self.name!=other.name
	
	def refresh_morgan(self):
		""" 
		Update the atom's morgan index
		:param self: this atom 
		"""
		self.morgan = self.new_morgan
		self.new_morgan=0
	
	def include_morgan_in_name(self):
		""" 
		Include the morgan index into the atom's name
		:param self: this atom
		"""
		self.name=self.old_name+str(self.morgan)
