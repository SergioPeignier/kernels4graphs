class bond:
	def __init__(self,covalence):
		""" creat a bond object. Covalence holds for the covalence of the bond, i.e., 
		the number of pairs of electron  shared by the atoms involved in the bond. 
		:param self: bond
		:param covalence : bond's covalence (it is also the bond's label)
		"""
		self.covalence=covalence
	
	def __str__(self):
		""" Overwrite the print function
		:param self: bond object
		"""	
		return "bond: "+str(self.covalence)	
	
	def __eq__(self,other):
		""" Overwrite the equality function
		:param self: bond
		:param other: another bond
		:type other: Bond obect 
		"""
		return self.covalence==other.covalence
	
	def __ne__(self,other):
		""" Overwrite the inequality function
		:param self: bond
		:param self: other bond
		:type other: Bond Object
		"""
		return self.covalence!=other.covalence