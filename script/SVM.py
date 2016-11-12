import numpy as np
import math
import cvxopt
import cvxopt.solvers as cvxopt_solvers

class SVM:
	def __init__(self,c,cutoff,activation):
		self.c=c
		self.cutoff=cutoff
		self.separator=[]
		self.activation=activation
	
	def answer(self, K_x):
		y_pred=0
		for i,K_x_i in enumerate(K_x) :
			y_pred+=self.separator[i]*K_x_i
		y_pred#+=self.b###################################################
		return self.activation(y_pred)
		
	def predict(self, K_pred_in):
		Y_pred=[]
		if type(Y_pred) != type([]):
			K_pred = K_pred_in.tolist()
		else:
			K_pred = K_pred_in
		for K_x in K_pred:
			Y_pred.append(self.answer(K_x))
		return np.array(Y_pred)	

	def to_sv(self,ranks,y):
		i=0
		sv = ranks.tolist()
		while i < len(sv):
			if sv[i]<self.cutoff:
				sv[i] = 0
			sv[i] *= y[i] 
			i+=1
		self.separator=sv
		return sv

	def fit(self,K_learn,y):
		self.K_matrix = K_learn
		n_samples=len(self.K_matrix)
		
		P=np.outer(y,y) * np.array(self.K_matrix)
		
		P = cvxopt.matrix(P,tc='d')
	
		q = cvxopt.matrix(np.ones(n_samples) * -1,tc='d')
		
		A = cvxopt.matrix(y, (1,n_samples),tc='d')
		
		b = cvxopt.matrix(0.0,tc='d')
	
		G = cvxopt.matrix(np.concatenate((np.diag(np.ones(n_samples))*-1,np.diag(np.ones(n_samples)))),tc='d')
		
		h = cvxopt.matrix(np.concatenate((np.zeros(n_samples),np.ones(n_samples)*self.c)),tc='d')
		
		# Solve QP problem
		solution = cvxopt_solvers.qp(P, q, G, h, A, b)
 
		# Lagrange multipliers
		self.to_sv(np.ravel(solution['x']),y)
		
		return self.separator
	
