import numpy as np
import networkx as nx
from kernel import kernel as kernel

class kernel_substrees_walk(kernel):
	def __init__(self,n=6):
		kernel.__init__(self)
		self.n=n
		
	def set_n(self,n):
		self.n = n	
	
	def cat_edge_node(self,node_id,node_father_id):
		bond=self.mol1.edges[self.mol1.edge_to_key([node_id,node_father_id])]
		return str(bond.covalence)+"-"+self.mol1.nodes[node_id].name
	
	def dico_of_neighbors(self,me,neighbors):
		dico={}
		for n in neighbors:
			key = self.cat_edge_node(n[0],me[0])
			if key not in dico.keys():
				dico[key] = []
			dico[key].append(n)
		return dico
		
	def not_visited(self,to_check,visited):
		not_visited=[]
		for v in to_check:
			if v not in visited:
				not_visited.append(v)
		return not_visited
		
	def tau_n(self,v,n,visited):
		
		neighboors=self.mol1_x_mol2.neighbors(v)
		if n > 1:
			usefull_neighboors=self.not_visited(neighboors,visited)
		else:return 1
		visited.append(v)
		visited2=visited[:]
		if usefull_neighboors:
			tau_tot = 0
			usefull_neighboors_dico = self.dico_of_neighbors(v,usefull_neighboors)
			for key in usefull_neighboors_dico.keys():
				tau_local=1
				for v_i in usefull_neighboors_dico[key]:
					visited2 = visited[:]
					p_t_mol1 = self.mol1.P_t(v[0], v_i[0])
					p_t_mol2 = self.mol2.P_t(v[1], v_i[1])
					lambd = p_t_mol1*p_t_mol2
					tau_child=lambd*self.tau_n(v_i,n-1,visited2)
					if tau_child:
						tau_local *= tau_child
					else:
						tau_local = 0
						break 
				tau_tot += tau_local
			return tau_tot
		else:return 0
		
	#compute the kernel function	
	def K(self,mol1,mol2):
		self.load(mol1,mol2)
		self.mol1_x_mol2=self.Mol_X_Mol(self.mol1,self.mol2)
		
		self.k=0
		size=self.n
		for atomes in self.mol1_x_mol2.nodes():
			tau=self.mol1.P_s(atomes[0])*self.mol2.P_s(atomes[1])
			visited = []
			self.k += tau*self.tau_n(atomes,size,visited)
		visited=[]
		return self.k
'''
import molecule as mc
mol1=["C","C","N","N"]
ed1=[(0,1),(1,2),(1,3)]
co1=[1,1,1]
mol2=["C","C","N","C","N"]
ed2=[(0,1),(1,2),(2,3),(2,4),(3,4)]
co2=[1,1,1,1,1]
comp1=mc.molecule(mol1,ed1,co1,0.1)
comp2=mc.molecule(mol2,ed2,co2,0.1)
k=kernel_substrees_walk(5)
print k.K(comp1,comp2)
'''



