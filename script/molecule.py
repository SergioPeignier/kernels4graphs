import networkx as nx
import atom
import bond
class molecule:
	def __init__(self, nodes, edges, covalances, p_q=0.1):#0.1
		self.G=nx.Graph()#Graph of molecules
		self.nodes=[]#dictionary of atoms (nodes)
		self.edges={}
		self.nodes=self.nodes_to_atoms(nodes)
		self.nodes_and_edges_2_graph(nodes, edges)
		self.set_deg_atoms()
		self.edges_to_bonds(edges,covalances)
		
		self.p_q=p_q#probability param => no transition probability
		self.p_o=1./(len(nodes))#initial probability distribution :> uniform distribution
		self.p_t_coef=(1.-p_q)#*1./p_q#transition probability coeficient 
		self.p_s=self.p_o * self.p_q#emission probability
	
	def nodes_to_atoms(self,nodes):
		"""Creates an atom array out of the atoms names
		:param nodes: names of atoms
		:type nodes: list of strings
		"""
		for node in nodes:
			self.nodes.append(atom.atom(node))
		return self.nodes
		
	def edge_to_key(self,edge):
		"""Converts pair of nodes (edge) into an identifier for a dictionary
		:param edge: List of atoms names forming an edge
		:type edge: List of Strings	
		"""
		edge=[edge[0],edge[1]]
		edge.sort()
		return str(edge[0])+"-"+str(edge[1])
		
	def edges_to_bonds(self,edges,covalances):
		""" Converts list of edges into bong objects
		:param edges: list of pairs of nodes names
		:param covalances: list of covalance values for each link
		:type edge: List
		:type covalence: List
		"""
		for i,edge in enumerate(edges):
			key=self.edge_to_key(edge)
			self.edges[key]=bond.bond(covalances[i])
		
	def nodes_and_edges_2_graph(self,nodes,edges):
		"""Makes a graph representing the molecule out of the nodes and edges list
		:param nodes: List of nodes
		:param edges: List of edges
		"""
		for i,node in enumerate(nodes):
			self.G.add_node(i,atom=node)
		self.G.add_edges_from(edges)
	
	def set_deg_atoms(self):
		"""
		Set the degrees for all the atoms of the molecules
		"""
		degrees=self.G.degree()
		for i,atom in enumerate(self.nodes):
			atom.deg=degrees[i]
	
	#set morgan indexes, T is the number of iterations
	def set_morgan_index(self,T):
		"""
		Set morgan index to each atom of the molecule (T is the number of iterations to update the morgan index)
		:param T: number of iterations to update the morgan index
		:type T: int
		"""
		i=0
		while i < T:
			for j, atom in enumerate(self.nodes):
				neighbors=self.G.neighbors(j)
				for neighbor in neighbors:
					atom.new_morgan+=self.nodes[neighbor].morgan
			for atom in self.nodes:
				atom.refresh_morgan()
			i+=1
		for atom in self.nodes:
				atom.include_morgan_in_name()
	
		
	def P_g(self,path):
		"""
		Compte the probability of a path given the probabilities received as attributes of the class
		:param path: path in the graph molecule
		:type path: list of edges
		"""
		ans=self.p_s
		i=1
		#print ans
		while i < len(path):
			#print self.P_t(path[i-1], path[i])
			ans*= self.P_t(path[i-1], path[i])
			i+=1
		return ans
	
	def P_s(self,node):
		"""
		Compute the emission probability (not used in the Kernel)
		:param node:
		:type node:
		"""
		return self.p_o * self.p_q

	def P_a(self,node):
		"""
		Compute the transition matrix probability
		:param node: node
		:type node: node id 
		"""
		if self.nodes[node].deg==0:	return 0
		else: return 1./self.nodes[node].deg
	
	def P_t(self, past, present):
		"""
		compute the transition probability
		:param past: previous node
		:param present: next node
		"""
		return self.p_t_coef*self.P_a(past)#*self.p_q
	
	def __str__(self):
		"""
		Overwrite print function
		"""
		txt="----------------------------------------\n"
		txt+= "---CHEMICAL COMPOUND DESCRIPTION---\n"
		txt+= "---nodes---\n"
		txt+= str(self.G.nodes()) + "\n"
		txt+= "--edges---\n"
		txt+= str(self.G.edges())+"\n"
		txt+= "---atoms details---\n"
		for i,atom in enumerate(self.nodes) :
			txt+= str(i)+"--> "+atom.__str__() +"\n"
		txt+="----------------------------------------\n"
		return txt