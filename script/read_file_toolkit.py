def usefull(w):
	return w!="" and w!="\n"
	
def read_line(line):
	if "\n" in line :line=line[0:len(line)-1]
	line=(" ".join(line.split("\t"))).split(" ")
	new_line=[]
	for w in line:
		if usefull(w):
			new_line.append(w)
	return new_line
		
def creat_local_txt(f,breaker="$"):
	txt=""
	while 1:
		a=f.readline()
		if "" == a:return -1
		if "$" == a[0]: break
		else: txt+=a
	return txt

def read_matrix(f,jump,to_num,stop):
	matrix = []
	i = 0
	while i < jump:
		f.readline()
		i+=1
	j=0	
	if to_num:
		while 1:
			line = read_line(f.readline())
			if not line or line[0] == stop:break 
			else: 
				matrix.append([])
				for a in line:
					matrix[j].append(float(a))
				j+=1
	else:
		while 1:
			line = read_line(f.readline())
			if not line or line[0] == stop:break 
			else: 
				matrix.append(line)
	return matrix
	
def read_row(f,start,to_num,stop):
	row = []
	read = False
	if to_num:
		while 1:
			line = read_line(f.readline())
			if not line or line[0] == stop:break
			if line[0] == start: read = True 
			else: 
				if read:
					for a in line:
						row.append(float(a))
				
	else:
		while 1:
			line = read_line(f.readline())
			if not line or line[0] == stop:break
			if line[0] == start: read = True 
			else: 
				if read:
					row = line
	return row
	
#f = open("./res/matrices/matrix_gram_bad_data_set_geometric_param_0.1","r")
#print(read_matrix(f,4,1,"***"))
#print(read_row(f,"***",1,""))
