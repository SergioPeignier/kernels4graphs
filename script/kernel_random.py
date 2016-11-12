import numpy as np
import networkx as nx
from kernel import kernel as kernel
class kernel_random(kernel):

	def __init__(self,p_q=0.1):
		kernel.__init__(self)
		self.p_q=p_q
		self.name="random_walks"
	def set_p_q(self,p_q):
		self.p_q = p_q
		
	def update_p_q(self):
		if self.mol1 != [] and self.mol2 != []:
			self.mol1.p_q=self.p_q
			self.mol2.p_q=self.p_q
		
	#generate the emissions vector	
	def PI_s(self):
		self.pi_s=[]
		for atomes in self.mol1_x_mol2.nodes():
			self.pi_s.append(self.mol1.P_s(atomes[0])*self.mol2.P_s(atomes[1]))
		self.pi_s=np.matrix(self.pi_s)
		return self.pi_s
		
	#generate the transition matrix	
	def PI_t(self):
		if self.mol1_x_mol2:
			self.pi_t=nx.adjacency_matrix(self.mol1_x_mol2)
		x_nodes=self.mol1_x_mol2.nodes()
		i=0
		while i < len(x_nodes):
			j=0
			while j < len(x_nodes):
				if self.pi_t[i,j] != 0:
					p_t_mol1=self.mol1.P_t(x_nodes[i][0], x_nodes[j][0])
					p_t_mol2=self.mol2.P_t(x_nodes[i][1], x_nodes[j][1])
					self.pi_t[i,j]=p_t_mol1*p_t_mol2
				j+=1
			i+=1
		return self.pi_t
		
	#compute the kernel function	
	def K(self,mol1,mol2):
		
		self.load(mol1,mol2)
		self.update_p_q()
		self.update_p_q()
		self.mol1_x_mol2=self.Mol_X_Mol(self.mol1,self.mol2)
		if self.len_prod_graph==0:return 0
		
		self.pi_s=self.PI_s()
		self.pi_t=self.PI_t()	
		
		I=np.ones(self.len_prod_graph)
		
		pi_t_power=np.eye(self.len_prod_graph)-self.pi_t
		pi_t_power=pi_t_power.I
	
		self.k=np.dot(np.dot(self.pi_s,pi_t_power),I)
		
		self.k=self.k[0,0]
	
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
k=kernel_random(0.1)
print k.K(comp1,comp2)
'''
