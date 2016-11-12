from read_file_toolkit import *
from stat_toolkit import *
import numpy as np
import SVM
import kernel_matrix
import kernel_random
import kernel_geometric
import kernel_nth
import kernel_subtrees
import kernel_subtrees_walk
import os
import sklearn
from sklearn import svm


name_ex = "./models/matrices/matrix_gram_bad_data_set_geometric_param_0.1"
params = {"radial":0.1,"normalized":0,"centered":0,"c":100,"cutoff":0.1}

def main_from_matrix(name, file_to_save, params):
	f = open(name,"r")
	M = read_matrix(f,4,1,"***")
	f.close()
	f = open(name,"r")
	Y = read_row(f,"***",1,"")
	f.close()

	k = ""
	kern = kernel_matrix.kernel_matrix(k,np.sign,params["radial"])
	#clf = svm.SVC(kernel = "precomputed",C=params["c"])
	clf = SVM.SVM(params["c"],params["cutoff"],np.sign)

	res = cross_validation_matrix(clf, kern, M, 10, params, Y)
	txt = statistics_to_txt(res,0,"roc_x roc_y score")
	params_learing_str = "_results_using_radial:"+str(params["radial"])+"_normalized:"+str(params["normalized"])+"_centered:"+str(params["centered"])+"_c:"+str(params["c"])+"_cutoff:"+str(params["cutoff"])
	g = open(file_to_save+"/"+name.split("/").pop()+params_learing_str,"w")
	g.write(txt)

matrices = os.listdir("./models/matrices/max")
print len(matrices)
for m in matrices:
	#try:
		#if "subtree" in m and "0.6" not in m and "." in m:
	print m
	name_dossier = "./results/results/max/2gaussian"
	params["radial"] = 0.0001
	params["normalized"] = 0
	params["centered"] = 0
	main_from_matrix("./models/matrices/max/"+m, name_dossier , params)
	#except:
	#	print "problem"
	#	pass	

'''
matrices = os.listdir("./res/matrices/max")

for m in matrices:
	
	name_dossier = "./res/results/max"
	params["radial"] = 3
	params["normalized"] = 0
	params["centered"] = 0
	main_from_matrix("./res/matrices/"+m, name_dossier , params)
	
'''
"""
c = 1
while c < 201:
	sigma = 0.002
	while sigma < 201:
		try:
			params["radial"] = sigma
			params["c"] = c
			main_from_matrix(name_ex , "./res/results/influence_of_sigma_c" , params)
			
		except:
			pass
		sigma += 2
	c+=10
"""
'''
c = 1
while c < 1001:
	try:
		params["radial"] = 0
		params["normalized"] = 0
		params["centered"] = 0
		params["c"] = c
		main_from_matrix(name_ex , "./res/results/influence_of_c" , params)
		
	except:
		pass
	c+=10
'''

#main_from_matrix(name_ex,"./res/results",params)
