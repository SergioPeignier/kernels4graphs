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

#loading data
cancer=setof.set_of_molecules()
nocancer=setof.set_of_molecules()
cancer.molecules_to_set("../data/toxic/cancer_nocancer_large/cancer.sdf",-1)
print "cancer loaded"
nocancer.molecules_to_set("../data/toxic/cancer_nocancer_large/nocancer.sdf",1)	
print "no cancer loaded"
print "seting to morgan"
cancer.set_morgan_to_all_molecules(1)
print len(cancer.my_set)
print len(nocancer.my_set)
nocancer.set_morgan_to_all_molecules(1)
print "data ok"
learn_set=setof.set_of_molecules()
test_set=setof.set_of_molecules()
learn_set.my_set = cancer.my_set[0:1000] + nocancer.my_set[0:1000]

random.shuffle(learn_set.my_set)
#test_set.my_set = cancer.my_set[100:150] + nocancer.my_set[300:350]
test_set.my_set = cancer.my_set[1000:1200] + nocancer.my_set[1000:1200]
random.shuffle(test_set.my_set)


#training and predicting

Y=tools_general.extract_column(learn_set.my_set,1)
kern = kernel_random.kernel_random(0.1)
#clf = svm.SVC(kernel = "precomputed")
param_kernel = {"morgan":1}
params = {"radial":0,"normalized":0,"centered":0,"c":550,"cutoff":0.1}#radial 0.002 ok
clf = SVM.SVM(params["c"],params["cutoff"],np.sign)
kern_perc=kernel_matrix.kernel_matrix(kern,np.sign,params["radial"])
gram = stat_toolkit.training(clf,kern_perc,learn_set,params)

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
