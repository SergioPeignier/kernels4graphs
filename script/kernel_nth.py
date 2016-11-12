import numpy as np
import networkx as nx
from kernel import kernel as kernel
class kernel_nth(kernel):
	def __init__(self,n=10):
		kernel.__init__(self)
		self.n=n
		self.name="power_adjacency"
	def set_n(self,n):
		self.n = n
		
	def K(self,mol1,mol2):
		self.load(mol1,mol2)
		self.mol1_x_mol2=self.Mol_X_Mol(self.mol1,self.mol2)
		A=nx.adjacency_matrix(self.mol1_x_mol2)
		A=np.matrix(A)
		size=len(A)
		An=np.linalg.matrix_power(A,self.n)
		I=np.matrix(np.ones(size))
		It=np.transpose(I)
		self.k=np.dot(np.dot(I,An),It)
		return self.k[0,0]
'''
import molecule as mc
mol1=["C","C","N","N"]
ed1=[(0,1),(1,2),(1,3)]
co1=[1,1,1]
mol2=["C","C","N","C","N"]
ed2=[(0,1),(1,2),(2,3),(2,4)]
co2=[1,1,1,1]
comp1=mc.molecule(mol1,ed1,co1,0.1)
comp2=mc.molecule(mol2,ed2,co2,0.1)
k=kernel_nth(3)
print k.K(comp1,comp2)
'''
