import numpy as np
import networkx as nx

class kernel:
	def __init__(self):
		self.mol1=[]
		self.mol2=[]
		self.mol1_x_mol2=0
		self.k=0
		self.pi_s=[]
		self.pi_t=[]
		self.len_prod_graph=0
		self.morgan = 0

	def set_morgan(self,morgan):
		self.morgan = morgan
		
	def  load(self,mol1,mol2):
		self.mol1=mol1
		self.mol2=mol2

	def change_morgan(self,morgan):
		if self.morgan != morgan:
			if self.mol1 != [] and self.mol2 != []:
				self.morgan = morgan
				self.mol1.set_morgan_to_all_molecules(morgan)
				self.mol2.set_morgan_to_all_molecules(morgan)
				
	#compute the graph product
	def Mol_X_Mol(self,mol1,mol2):	
		self.mol1_x_mol2 = nx.algorithms.product.tensor_product(mol1.G,mol2.G)
		for node in self.mol1_x_mol2.nodes():
			if self.mol1.nodes[node[0]] != self.mol2.nodes[node[1]]: self.mol1_x_mol2.remove_node(node)
		
		for edge in self.mol1_x_mol2.edges():
			bond1=self.mol1.edges[self.mol1.edge_to_key([edge[0][0],edge[1][0]])]
			bond2=self.mol2.edges[self.mol2.edge_to_key([edge[0][1],edge[1][1]])]
			if bond1 != bond2: self.mol1_x_mol2.remove_edge(edge[0],edge[1])
		self.len_prod_graph=len(self.mol1_x_mol2.nodes())
		
		return self.mol1_x_mol2
