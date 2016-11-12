import numpy as np
import math

class kernel_matrix:
	#training=[[molecule,label]...]
	def __init__(self,kernel,activation,sigma=2,c=0.1):
		self.kernel=kernel
		self.K_matrix=[]
		self.training=[]
		self.activation=activation
		self.sigma=sigma
		
	def K_computation(self,training):
		self.training=training
		n_samples=len(training)
		self.K_matrix = np.zeros((n_samples, n_samples))
		i=0
		while i < n_samples:
			j=0
			self.K_matrix[i,i] = self.kernel.K(training[i][0], training[i][0]) 
			while j < i:
				self.K_matrix[i,j] = self.kernel.K(training[i][0], training[j][0])
				self.K_matrix[j,i] = self.K_matrix[i,j]
				j+=1
			print i
			i+=1
		return self.K_matrix
	
	def normalize_k(self):
		normalized=[]
		for i,x1 in enumerate(self.K_matrix):
			normalized.append([])
			for j,x2 in enumerate(self.K_matrix):
				normalized[i].append(self.K_matrix[i,j]*1./math.sqrt(self.K_matrix[i,i]*self.K_matrix[j,j]+0.00000000001))
		normalized=np.matrix(normalized)
		return normalized
		
	def center_k(self):
		n = len(self.K_matrix)
		Z=(np.matrix(np.identity(n))-1./n*np.matrix(np.ones((n,n))))
		Kc=np.dot(np.dot(Z,self.K_matrix),Z)
		return Kc

	def K_matrix_to_radial_bases(self):
		K_radial=[]
		for i,x1 in enumerate(self.K_matrix):
			K_radial.append([])
			for j,x2 in enumerate(self.K_matrix):
				K_radial[i].append(self.k_to_radial(self.K_matrix[i,i],self.K_matrix[i,j],self.K_matrix[j,j]))
		K_radial=np.matrix(K_radial)
		return K_radial

	def D2(self,kii,kij,kjj):
		return kii - 2*kij + kjj
		
	def k_to_radial(self,kii,kij,kjj):
		return math.exp(-self.D2(kii,kij,kjj)*1./(2*self.sigma))
		
	def k_func_radial(self,x,x_base,i):
		kbase=self.K_matrix[i,i]
		kx=self.kernel.K(x[0],x[0])
		k=self.kernel.K(x_base[0],x[0])
		return self.k_to_radial(kbase,k,kx)
	
	def k_func_norm(self,x,x_base,i):
		kbase=self.K_matrix[i,i]
		kx=self.kernel.K(x[0],x[0])
		k=self.kernel.K(x_base[0],x[0])
		return k*1./math.sqrt(kx*kbase)
	
	def k_func_simple(self,x,x_base,i):
		return self.kernel.K(x_base[0],x[0])
		
	def k_answer(self,x,x_base_i,func):
		return func(x,x_base,i)
		
	def k_answer_many(self,test_data,func,sv):
		K_ans=[]
		if sv == "x":
			for i,indiv in enumerate(test_data):
				K_ans.append([])
				for j, k in enumerate(self.training):
					K_ans[i].append(func(indiv,self.training[j],j))
		else:
			for i,indiv in enumerate(test_data):
				K_ans.append([])
				for j, k in enumerate(self.training):
					if sv[j]:
						K_ans[i].append(func(indiv,self.training[j],j))
					else:
						K_ans[i].append(0)
					
		K_ans=np.matrix(K_ans)
		return K_ans
		
	def diversity_mesure(self,training):
		if self.K_matrix==[]:self.K_matrix=self.K_computation(training)
		D=0
		i=0
		local=0
		while i < len(self.K_matrix):
			D+=self.K_matrix[i,i]
			j=0
			while j < i:
				local+=2*self.K_matrix[i,j]
				j+=1
			i+=1
		local*=1./(len(self.K_matrix)**2)
		D*=(1./len(self.K_matrix)-1./len(self.K_matrix)**2)
		D-=local
		return D
	
	
