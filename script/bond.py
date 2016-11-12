class bond:
	#covalence=1,2 ...
	#further energy, length ....
	""" creat a bond object
		:param self: bond
		:param covalence : bond's covalence (it is also the bond's label)
		"""
	def __init__(self,covalence):
		self.covalence=covalence
	""" surcharge the print function
		:param self: bond
		"""	
	def __str__(self):
		return "bond: "+str(self.covalence)	
	""" surcharge the equality function
		:param self: bond
		:param other: another bond 
		"""
	def __eq__(self,other):
		return self.covalence==other.covalence
	""" surcharge the inequality function
		:param self: bond
		:param self: other bond
		"""
	def __ne__(self,other):
		return self.covalence!=other.covalence
'''		
b=bond(1)
c=bond(1)
d=bond(2)
print b==c
print b==d
print b!=d
print b!=c
'''
