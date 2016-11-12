import kernel_matrix
import numpy as np
import random
def statistics_to_txt(cross,column=0,legend=""):
	txt=""
	if column:
		for a in zip(*cross):
			for value in a:
				txt+=str(value)+" "
			txt+="\n"
	else:
		for a in cross:
			for value in a:
				txt+=str(value)+" "
			txt+="\n"
	return legend+"\n"+txt
	
def splitter(data,i,parts,which="not random"):
	ans={"test":[],"training":[]}
	if which == "random":
		data.shuffle()
		cut=1*len(data)/parts
		ans["test"]=data[:cut]
		ans["training"]=data[cut:]
	else: 
		init=int(i*len(data)/parts)
		end=int((i+1)*len(data)/parts)
		if end>=len(data): 
			ans["test"]=data[init:]
			ans["training"]=data[0:init]
		else: 
			ans["test"]=data[init:end]
			ans["training"]=data[0:init]+data[end:]
	return ans

def shuffle_matrix(matrix , Y):
	new_order=range(len(matrix))
	random.shuffle(new_order)
	new_matrix=[]
	y = []
	for i in new_order:
		new_matrix.append(matrix[i])
		y.append(Y[i])
	answer_matrix = []
	for i,col in enumerate(new_matrix):
		answer_matrix.append([])
		for j in new_order:
			answer_matrix[i].append(col[j])
	return {"matrix":answer_matrix,"Y":y}
	
def extract_column(training,col):
	label=[]
	for x in training:
		label.append(x[col])
	return label

#[TP,FP,TN,FN]
def error(pred,ans):
	
	if ans==1:
		if pred == ans:return [1,0,0,0]
		else: return [0,0,0,1]
	if ans==-1:
		if pred == 1:return [0,1,0,0]
		else: return [0,0,1,0]
	if ans == 0:
		return 0
		
def ROC(prediction,ans):
	#print prediction
	#print ans
	t_p=0
	f_p=0
	t_n=0
	f_n=0
	for i,x in enumerate(ans):
		my_error=error(prediction[i],x)
		t_p+=my_error[0]
		f_p+=my_error[1]
		t_n+=my_error[2]
		f_n+=my_error[3]

	ans=[t_p*1./(t_p+f_n+0.1),t_n*1./(t_n+f_p+0.1),(t_p+t_n)*1./(f_p+t_n+t_p+f_n+0.1)]
	return ans
#clf = svm.SVC(kernel = "precomputed")
def compute_gram(kern_perc,params):
	
	if params["centered"]:
		return kern_perc.center_k()
	if params["normalized"]:
		return kern_perc.normalize_k()
	if params["radial"]:
		return kern_perc.K_matrix_to_radial_bases()
	return kern_perc.K_matrix
	
def training(clf,kern_perc,learn_set,params):
	kern_perc.K_computation(learn_set.my_set)
	gram = compute_gram(kern_perc,params)
	Y = extract_column(learn_set.my_set,1)
	clf.fit(gram,Y)
	return gram
	
def predicting(clf,kern_perc,test_set,params,sv):	
	if params["radial"]:
		print "using radial to predict"
		Kans = kern_perc.k_answer_many(test_set.my_set,kern_perc.k_func_radial,sv)
	else:
		if params["normalized"]:
			print "using normalized to predict"
			Kans = kern_perc.k_answer_many(test_set.my_set,kern_perc.k_func_norm,sv)
		else:
			print "using simple to predict"
			Kans = kern_perc.k_answer_many(test_set.my_set,kern_perc.k_func_simple,sv)
	prediction=clf.predict(Kans)
	return prediction			
			
def cross_validation(clf,kernel_func,training_set,which_splitter,nb,params):#p=answer function parameter=function to compute K(x,x')
	i=0
	error=0
	cross_val_error=[]
	while i < nb:
		#generate data sets
		splitted = splitter(training_set,i,nb,which_splitter)
		test = splitted["test"]
		train = splitted["training"]
		#creating the kernel
		kern_perc=kernel_matrix.kernel_matrix(kernel_func,np.sign,params["sigma"])
		#training the model
		training(clf,kernel_func,train,params)
		#predicting
		prediction = predicting(clf,kern_perc,test,params)
		#getting true answers
		Y_test = extract_column(test,1)
		#computing and saving scores
		cross_val_error.append(ROC(prediction,Y_test))
		
		i+=1
		print i
	return cross_val_error
	
def set_matrixes (i,data,parts,Y):
	init=int(i*len(data)/parts)
	end=int((i+1)*len(data)/parts)
	
	ans = {"test" : [], "training" : [] ,"Y_test" : [], "Y_training" : []}
	if end>=len(data): 
		i = init
		while i < len(data):
			ans["test"].append(data[i][0:init])
			ans["Y_test"].append(Y[i])
			i += 1
		i = 0
		while i < init:
			ans["training"].append(data[i][0:init])
			ans["Y_training"].append(Y[i])
			i += 1
	else: 
		
		i = init
		while i < end:
			ans["test"].append(data[i][0:init]+data[i][end:])
			ans["Y_test"].append(Y[i])
			i += 1
		i = 0
		while i < init:
			ans["training"].append(data[i][0:init]+data[i][end:])
			ans["Y_training"].append(Y[i])
			i += 1
		i = end
		while i < len(data):
			ans["training"].append(data[i][0:init]+data[i][end:])
			ans["Y_training"].append(Y[i])
			i += 1
	return ans		
	
def predicting_from_matrix(clf, matrix):					
	prediction=clf.predict(matrix)	
	return prediction
	
def training_form_matrix(clf, K, Y):
	clf.fit(K,Y)
	return K		
	
def cross_validation_matrix(clf, kern_perc, K_M, nb, params, Y):
	kern_perc.K_matrix = np.matrix(K_M)
	K = np.matrix(K_M)
	if params["radial"]:
		print "radializando"
		K = kern_perc.K_matrix_to_radial_bases()
	else:
		if params["normalized"]:
			print "normando"
			K = kern_perc.normalize_k()
		else:
			if params["centered"]:
				print "centrando"
				K = kern_perc.center_k()
	#data_set = shuffle_matrix(K_M, Y)
	matrix = K.tolist()#data_set["matrix"]
	y =  Y#data_set["Y"]
	cross_val_error=[]
	i = 0
	while i < nb:
		cutted = set_matrixes (i, matrix, nb, Y)
		training_form_matrix(clf, cutted["training"], cutted["Y_training"])
		prediction = predicting_from_matrix(clf, cutted["test"])
		cross_val_error.append(ROC(prediction,cutted["Y_test"]))
		i+=1
	mean=[0,0,0]
	for res in cross_val_error:
		for j in range(3):
			mean[j]+=res[j]
	for j in range(3):
		mean[j]*=1./len(cross_val_error)
	cross_val_error.append(mean)
	return cross_val_error

		
def dico_to_txt(txt_in,dico):
	txt = txt_in
	for k in dico.keys():
		txt+="|"+k+" : "+str(dico[k])+"|"
	return txt
		
def generate_legend_for_results(params_SVM,param_kernel,pred_or_test,more_details=""):
	svm_txt = dico_to_txt("The SVM parameters are: ",params_SVM)
	kern_txt = dico_to_txt("The kernel fixed parameters are: ",param_kernel)
	return pred_or_test+'\n'+svm_txt+"\n"+kern_txt+"\n"+more_details
	

def eval_param(clf,kernel_func,learn_set,test_set,params,set_param,where):
	Y=extract_column(learn_set.my_set,1)
	Ytest = extract_column(test_set.my_set,1)
	result = {"training":[],"prediction":[]}
	for val in where:
		try:
			set_param(val)
			kern_perc=kernel_matrix.kernel_matrix(kernel_func,np.sign,params["radial"])
			gram = training(clf,kern_perc,learn_set,params)
			
			pred = predicting(clf,kern_perc,test_set,params)
			training_pred=clf.predict(gram)
			
			result["training"].append([val]+ROC(training_pred,Y))
			result["prediction"].append([val]+ROC(pred,Ytest))
		except:
			pass
	return result
