def extract_column(training,col):
	label=[]
	for x in training:
		label.append(x[col])
	return label
