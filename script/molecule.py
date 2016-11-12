import networkx as nx
import atom
import bond
class molecule:
	#self.G=molecular graph NX
	#self.nodes=dictionary of atoms (nodes)
	def __init__(self, nodes, edges, covalances, p_q=0.1):#0.1
		self.G=nx.Graph()
		self.nodes=[]
		self.edges={}
		self.nodes=self.nodes_to_atoms(nodes)
		self.nodes_and_edges_2_graph(nodes, edges)
		self.set_deg_atoms()
		self.edges_to_bonds(edges,covalances)
		
		self.p_q=p_q#probability param :> no transition
		self.p_o=1./(len(nodes))#initial probability distribution :> uniform distribution
		self.p_t_coef=(1.-p_q)#*1./p_q#transition probability coeficient 
		
		self.p_s=self.p_o * self.p_q#emission probability
	
	
	
	#uses an atoms names array to creat an atoms array
	def nodes_to_atoms(self,nodes):
		for node in nodes:
			self.nodes.append(atom.atom(node))
		return self.nodes
		
	def edge_to_key(self,edge):
		edge=[edge[0],edge[1]]
		edge.sort()
		return str(edge[0])+"-"+str(edge[1])
		
	def edges_to_bonds(self,edges,covalances):
		for i,edge in enumerate(edges):
			key=self.edge_to_key(edge)
			self.edges[key]=bond.bond(covalances[i])
		
	#make a graph from nodes(atoms) and edges
	def nodes_and_edges_2_graph(self,nodes,edges):
		for i,node in enumerate(nodes):
			self.G.add_node(i,atom=node)
		self.G.add_edges_from(edges)
	
	#set all the molecules atoms degrees
	def set_deg_atoms(self):
		degrees=self.G.degree()
		for i,atom in enumerate(self.nodes):
			atom.deg=degrees[i]
	
	#set morgan indexes, T is the number of iterations
	def set_morgan_index(self,T):
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
	
		
	#compute a paths probability
	def P_g(self,path):
		ans=self.p_s
		i=1
		#print ans
		while i < len(path):
			#print self.P_t(path[i-1], path[i])
			ans*= self.P_t(path[i-1], path[i])
			i+=1
		return ans
	
	#compute the emission probability (useless here)
	def P_s(self,node):
		return self.p_o * self.p_q

	#compute the transition matrix probability
	def P_a(self,node):
		if self.nodes[node].deg==0:	return 0
		else: return 1./self.nodes[node].deg
	
	#compute the transition probability
	def P_t(self, past, present):
		return self.p_t_coef*self.P_a(past)#*self.p_q
	
	#print
	def __str__(self):
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

'''
mol1=["C","C","C","O"]
mol2=["C","C","C","O"]
mol3=["H","C","O","Cl"]
mol4=["C","C","C","C","C","C","C","C","C","C","N","O","O"]

edgesmol1=[(0,1),(1,2),(2,3),(3,0)]
edgesmol2=[(0,1),(1,2),(2,3),(3,0)]
edgesmol3=[(0,1),(1,2),(1,3)]
edgesmol4=[(0,1),(2,0),(1,3),(2,4),(5,3),(4,5),(6,4),(7,5),(8,6),(9,7),(9,8),(10,9),(11,10),(12,10)]

comp1=molecule(mol1,edgesmol1,0.1)
comp2=molecule(mol2,edgesmol2,0.1)
comp3=molecule(mol3,edgesmol3,0.1)
comp4=molecule(mol4,edgesmol4,0.1)

print comp1
print comp2
print comp3
print comp4



comp4.set_morgan_index(1)
print comp4
comp4.set_morgan_index(1)
print comp4


m=nx.adjacency_matrix(comp4.G)
print m[1,1]
x_nodes=comp4.G.nodes()
i=0
while i < len(x_nodes):
	j=0
	while j < len(x_nodes):
		if m[i,j] :
			m[i,j]=comp4.P_t(i,j)
		j+=1
	i+=1


print m
'''		
