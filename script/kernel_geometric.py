import numpy as np
import networkx as nx
from kernel import kernel as kernel

class kernel_geometric(kernel):
	def __init__(self,b=0.6):
		kernel.__init__(self) 
		self.b=b
		self.name="geometric"
	def set_b (self,b):
		self.b = b
		
	#generate the emissions vector	
	def PI_s(self):
		self.pi_s=[]
		for atomes in self.mol1_x_mol2.nodes():
			self.pi_s.append(1)
		self.pi_s=np.matrix(self.pi_s)
		return self.pi_s
		
	#generate the transition matrix	
	def PI_t(self):
		self.pi_t=nx.adjacency_matrix(self.mol1_x_mol2)
		i=0
		while i < self.len_prod_graph:
			j=0
			while j < self.len_prod_graph:
				if self.pi_t[i,j] != 0:
					self.pi_t[i,j]=self.b
				j+=1
			i+=1
		return self.pi_t
		
	def K(self,mol1,mol2):
		
		self.load(mol1,mol2)
		self.mol1_x_mol2=self.Mol_X_Mol(self.mol1,self.mol2)
		
		self.pi_s=self.PI_s()
		self.pi_t=self.PI_t()
		
		size=len(self.mol1_x_mol2.nodes())
		if size==0:return 0
		
		I=np.ones(size)
		pi_t_power=np.eye(size)-self.pi_t
		pi_t_power=pi_t_power.I
	
		self.k=np.dot(np.dot(self.pi_s,pi_t_power),I)
		self.k=self.k[0,0]
		if self.k<0:self.k=0# beware this kernel may not converge (by definition of this kernel)!
		return self.k
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
k=kernel_geometric(0.01)
print k.K(comp1,comp2)
'''
