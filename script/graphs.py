import networkx as nx
import numpy as np


class atom:
	#name='C'#carbon ...
	#deg=2#number of neigboors
	#morgan=3 morgan coeficient	
	def __init__(self,name):
		self.name=name
		self.deg=0
		self.morgan=1
		self.new_morgan=0
	
	def __str__(self):
		return "atom: "+self.name+" | deg: "+str(self.deg)+" | morgan: "+str(self.morgan)	
	def __eq__(self,other):
		return self.name==other.name
	def __ne__(self,other):
		return self.name!=other.name
	def refresh_morgan(self):
		self.morgan = self.new_morgan
		self.new_morgan=0
	def include_morgan_in_name(self):
		self.name=self.name+str(self.morgan)

class molecule:
	#self.G=molecular graph NX
	#self.nodes=dictionary of atoms (nodes)
	def __init__(self, nodes, edges, p_q):
		self.G=nx.Graph()
		self.nodes=nodes
		self.nodes_and_edges_2_graph(nodes, edges)
		self.set_deg_atoms()
		
		self.p_q=p_q#probability param :> no transition
		self.p_o=1./(len(nodes))#initial probability distribution :> uniform distribution
		self.p_t_coef=(1.-p_q)#*1./p_q#transition probability coeficient 
		
		self.p_s=self.p_o * self.p_q#emission probability
		
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
	'''
	#set morgan indexes, T is the number of iterations
	def set_morgan_index(self,T):
		i=0
		while i < T:
			for i, atom in enumerate(self.nodes):
				#neighboors=obtener vecinos (seguro existe metodo)
				for neighboor in neighboors:
					atom.new_morgan+=neighboors.morgan
			for atom in self.nodes:
				atom.refresh_morgan()
			i+=1
		for atom in self.nodes:
				atom.include_morgan_in_name()
	'''
		
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
		txt=""
		txt+= "---CHEMICAL COMPOUND DESCRIPTION---\n"
		txt+= "---nodes---\n"
		txt+= str(self.G.nodes()) + "\n"
		txt+= "--edges---\n"
		txt+= str(self.G.edges())+"\n"
		txt+= "---atoms details---\n"
		for i,atom in enumerate(self.G.nodes(2)) :
			txt+= str(i)+"--> "+atom[1]["atom"].__str__() +"\n"
		return txt
#for the kernel 
	#def __mull__(self,other):
	#	ans = molecule([],[])
	#	
	#	return ans

class kernel:
	#self.G=molecular graph NX
	#self.nodes=dictionary of atoms (nodes)
	def __init__(self):
		self.mol1=[]
		self.mol2=[]
		self.mol1_x_mol2=0
		self.k=0
		self.pi_s=[]
		self.pi_t=[]
	#multiply 2 chemical compounds graphs and choose in these graph only the nodes with 2 sames atoms inside 
	def  load(self,mol1,mol2):
		self.mol1=mol1
		self.mol2=mol2
		 
	def Mol_X_Mol(self,mol1,mol2):
		self.mol1_x_mol2 = nx.algorithms.product.cartesian_product(mol1.G,mol2.G)
		#print self.mol1_x_mol2.edges()
		#print self.mol1_x_mol2.nodes()
		for node in self.mol1_x_mol2.nodes():
			if self.mol1.nodes[node[0]] != self.mol2.nodes[node[1]]: self.mol1_x_mol2.remove_node(node)
		return self.mol1_x_mol2
		
	#generate the emissions vector	
	def PI_s(self):
		for atomes in self.mol1_x_mol2:
			self.pi_s.append(self.mol1.P_s(atomes[0])*self.mol2.P_s(atomes[1]))
		self.pi_s=np.matrix(self.pi_s)
		return self.pi_s
		
	#generate the transition matrix	
	def PI_t(self):
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
		
	#generate a vector of a repetition of numbers
	def identical(self,number,iterations):
		vec=[]
		i=0
		while i < iterations:
			vec.append(number)
			i+=1
		return vec
	
	#compute the kernel function	
	def K(self,mol1,mol2):
		self.load(mol1,mol2)
		self.mol1_x_mol2=self.Mol_X_Mol(self.mol1,self.mol2)
		if self.pi_s==[]:self.pi_s=self.PI_s()
		if self.pi_t==[]:self.pi_t=self.PI_t()
		
		size=len(self.mol1_x_mol2.nodes())
		I=self.identical(1.,size)
		
		pi_t_power=np.identity(size)-self.pi_t
		#print pi_t_power
		pi_t_power=pi_t_power.I
		#print pi_t_power
		self.k=np.dot(self.pi_s,pi_t_power)
		#print self.k
		self.k=np.dot(self.k,I)
		#print "!!!!!!!!!!"
		#print self.k
		self.k=self.k[0,0]
		return self.k
	

class SVM:
	 #training=[[molecule,label]...]
	def __init__(self,kernel):
		self.kernel=kernel
		self.K_matrix=[]
		self.training=[]
		self.separator=[]
		
	def load_training(self,training):
		self.training = training
		
	def init_separator(self,training):
		self.separator=self.identical(0,len(training))
		self.separator=np.matrix(self.separator)
		
	def K_computation(self,training):
		for i,x1 in enumerate(training):
			self.K_matrix.append([])
			for j,x2 in enumerate(training):
				self.K_matrix[i].append(self.kernel.K(x1[0],x2[0]))
		self.K_matrix=np.matrix(self.K_matrix)
		return self.K_matrix
		
	def identical(self,number,iterations):
		vec=[]
		i=0
		while i < iterations:
			vec.append(number)
			i+=1
		return vec
	
	def answer(self,separator,training,x):
		y_pred=0
		for i,x_t in enumerate(training) :
			y_pred+=separator[0,i]*self.kernel.K(x_t[0],x[0])
		if y_pred>0 : return 1
		else: return -1
		
	def answer_many(self,separator, training, X):
		Y_pred=[]
		#print "anssss"
		#print separator
		#print training
		#print X
		for x in X:
			Y_pred.append(self.answer(separator,training,x))
		return Y_pred
		
	def learn(self,training,T):
		self.training = training
		self.K_matrix=self.K_computation(training)
		self.init_separator(training)
		i=0
		j=0
		#print "OKKKSS"
		#print self.separator
		#print self.K_matrix
		
		while i<len(training) and j<T:
			ans=training[i][1]*np.dot(self.separator,np.transpose(self.K_matrix[i,]))
			ans=ans[0,0]
		#	print "training"
		#	print training[i][1]
		#	print np.dot(self.separator,np.transpose(self.K_matrix[i,]))
			if ans > 0: i+=1
			else: self.separator[0,i]+=training[i][1]
			j+=1
		#print j
		return  self.separator
mol1=["C","C","C","O"]
mol2=["C","C","C","O"]
mol3=["H","C","O","Cl"]

edgesmol1=[(0,1),(1,2),(2,3),(3,0)]
edgesmol2=[(0,1),(1,2),(2,3),(3,0)]
edgesmol3=[(0,1),(1,2),(1,3)]

nodes1=[]
nodes2=[]
nodes3=[]

for a in mol1:
	nodes1.append(atom(a))
for a in mol2:
	nodes2.append(atom(a))
for a in mol3:
	nodes3.append(atom(a))

comp1=molecule(nodes1,edgesmol1,0.1)
comp2=molecule(nodes2,edgesmol2,0.1)
comp3=molecule(nodes3,edgesmol3,0.1)

#print comp1
#print comp2
#print comp3
#print comp3.P_g([0,1,2])

kern=kernel()
#kern.load(comp1,comp2)
#h=kern.Mol_X_Mol(kern.mol1,kern.mol2)
#print kern.PI_s()
#print kern.PI_t()
print kern.K(comp1,comp2)
#import Graphics

#np.array([100., 100.])
#a = np.matrix('1 2; 3 4')
#x = np.array([[1,2,3], [4,5,6]])
