import stat_toolkit
import random
import numpy as np
import math
import sklearn
from sklearn import svm
import kernel_matrix
import kernel_random
import kernel_geometric
import kernel_nth
import kernel_subtrees
import kernel_subtrees_walk
import set_of_molecules as setof
import SVM
import tools_general



random.seed(1)	

data=setof.set_of_molecules()
data.smiles_files_to_molecules("../data/cpdb_mutagens/mutagen.tab",1,2)
data.set_morgan_to_all_molecules(1)
random.shuffle(data.my_set)
learn_set=setof.set_of_molecules()
test_set=setof.set_of_molecules()
learn_set.my_set=data.my_set
#test_set.my_set=data.my_set[4:6]#len(data.my_set)-1]

Y=tools_general.extract_column(learn_set.my_set,1)

param_kernel = {"morgan":1}
params = {"radial":0,"normalized":0,"centered":0,"c":150,"cutoff":0.1}#radial 0.002 ok   #para otro c = 150 cool
pq=[0.7,0.8]
n=[2,4]
b=[0.2,0.3,0.4,0.5,0.6,0.7,0.8]
nst=[1,2,3,4,5,6,7,8,9]
lambs=[0.01,0.05,0.075,0.1,0.125,0.15,0.2,0.3]
lambs=[0.075,0.125]
#lambs = [0.25,0.35,0.45,0.7,0.8]
morgan = [3]
for p in nst:#lambs:
	try:
		kern = kernel_random.kernel_random(p)
		learn_set.set_morgan_to_all_molecules(p)
		#kern = kernel_nth.kernel_nth(p)
		#kern = kernel_geometric.kernel_geometric(p)
		#kern = kernel_subtrees.kernel_subtrees(p,0.1)#2,p)
		clf = SVM.SVM(params["c"],params["cutoff"],np.sign)
		kern_perc=kernel_matrix.kernel_matrix(kern,np.sign,params["radial"])
		gram = kern_perc.K_computation(learn_set.my_set)#stat_toolkit.training(clf,kern_perc,learn_set,params)
		header_test=stat_toolkit.generate_legend_for_results(params,param_kernel,"matrix","parameter teste: morgan = "+str(p))
		grama=gram.tolist()
		#print grama
		gram_txt = stat_toolkit.statistics_to_txt(grama,0,header_test)
		gram_txt += "***\n"
		for yi in Y:
			gram_txt += str(yi)+" "
		f = open("./models/matrices/matrix_gram_bad_data_set_"+kern.name+"_param_"+str(p),"w")
		f.write(gram_txt)
	except:
		print "pb"
		pass
'''
try:
	sv = clf.separator
except:
	sv = "x"

pred = stat_toolkit.predicting(clf,kern_perc,test_set,params,sv)
training_pred=clf.predict(gram)
Ytest = tools_general.extract_column(test_set.my_set,1)

print "error in prediction"
print stat_toolkit.ROC(pred,Ytest)
print "error in training"
print stat_toolkit.ROC(training_pred,Y)

'''
