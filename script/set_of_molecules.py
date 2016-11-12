import molecule
import read_file_toolkit as rf
from rdkit import Chem
from rdkit.Chem import AllChem



class set_of_molecules:
	def __init__(self,p_q=0.1):
		self.my_set=[]
		self.type_bonds={'SINGLE': 1,'DOUBLE': 2,'TRIPLE': 3, 'AROMATIC': 4}
		self.p_q=p_q

	def mol_to_nodes_and_edges(self, txt, class_label="", positive_class="active", negative_class="inactive"):
		txt=txt.split("\n")
		#print txt
		nodes=[]
		edges=[]
		covalance=[]
		head=rf.read_line(txt[3])
		n_nodes=int(head[0])
		n_edges=int(head[1])
		if n_nodes==0 or n_edges==0: return {"nodes":nodes,"edges":edges,"class":"","covalences":covalance}
		i=1
		while i<=n_nodes:
			line=rf.read_line(txt[3+i])
			#if len(line)<4: 
				#print line
				#print head
			nodes.append(line[3])
			i+=1
		i=1
		while i<=n_edges:
			line=rf.read_line(txt[3+n_nodes+i])
			edges.append((int(line[0])-1,int(line[1])-1))
			covalance.append(int(line[3]))
			i+=1
		i=3+n_nodes+n_edges
		if class_label=="":return {"nodes":nodes,"edges":edges,"class":"","covalences":covalance}
		
		while i<len(txt):
			#print class_label
			if  class_label in txt[i]:
				#print positive_class
				#print txt[i+1]
				#print negative_class
				if positive_class == txt[i+1][0:len(positive_class)]: return {"nodes":nodes,"edges":edges,"class":1,"covalences":covalance}
				if negative_class == txt[i+1][0:len(negative_class)] : return {"nodes":nodes,"edges":edges,"class":-1,"covalences":covalance}

			i+=1
		return {"nodes":nodes,"edges":edges,"class":"?","covalences":covalance}

	#for SMILES
	def smile_to_nodes_and_edges(self,smile):
		smile=smile[1:len(smile)-1]
		m = Chem.MolFromSmiles(smile,0)
		atoms = m.GetAtoms()
		bonds = m.GetBonds()
		nodes=[]
		edges=[]
		covalences=[]
		for a in atoms:
			nodes.append(a.GetSymbol())
		for b in bonds:
			covalences.append(self.type_bonds[b.GetBondType().__str__()])
			edges.append((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))
		return {"nodes":nodes,"edges":edges,"covalences":covalences}
		#[molecule.molecule(mol_loc["nodes"],mol_loc["edges"],mol_loc["covalences"]),mol_loc["class"]]
	
	#file in tab format
	def smiles_files_to_molecules(self,file_name,smile_pos,class_pos=None):
		compounds=[]
		f=open(file_name,"r")
		while 1:
			txt=rf.read_line(f.readline())
			if txt==[]:break
			mol_loc=self.smile_to_nodes_and_edges(txt[smile_pos]) 
			if mol_loc["nodes"]!=[] and mol_loc["edges"]!=[]:
				clas=txt[class_pos]
				if clas:
					clas=int(clas)
					if clas==0:clas=-1 
					compounds.append([molecule.molecule(mol_loc["nodes"],mol_loc["edges"],mol_loc["covalences"],self.p_q),clas])
				else:compounds.append([molecule.molecule(mol_loc["nodes"],mol_loc["edges"],mol_loc["covalences"],self.p_q),'?'])
		self.my_set=compounds
		return compounds
		
	#for MOL files 	
	def all_file_to_molecules(self, file_name,class_label="",pos_class="active",neg_class="inactive"):
		compounds=[]
		f=open(file_name,"r")
		while 1:
			txt=rf.creat_local_txt(f)
			if txt == -1: return compounds
			else:
				mol_loc=self.mol_to_nodes_and_edges(txt,class_label,pos_class) 
				#print mol_loc
				if mol_loc["nodes"]!=[] and mol_loc["edges"]!=[]:
					compounds.append([molecule.molecule(mol_loc["nodes"],mol_loc["edges"],mol_loc["covalences"],self.p_q),mol_loc["class"]])
		return compounds
	
	def all_files_to_molecules(self, files_names,class_label="",pos_class="active",neg_class="inactive"):
		compounds=[]
		for file_name in files_names:
			compounds+=all_file_to_molecules(file_name, class_label,pos_class,neg_class)
		return compounds
	
	def files_to_nodes_general(self, files_names,class_label="",pos_class="active",neg_class="inactive"):
		if type(files_names)==type(""): return self.all_file_to_molecules(files_names,class_label,pos_class,neg_class)
		if type(files_names)==type([]): return  self.all_files_to_molecules(files_names,class_label,pos_class,neg_class)
		
	
	def molecules_to_set(self, files_names,y):
		molecules = self.files_to_nodes_general(files_names)
		i=0
		while i < len(molecules):
			molecules[i][1]=y
			i+=1
		self.my_set=molecules
		return self.my_set
		
	def set_morgan_to_all_molecules(self,T=2):
		for i,molecule in enumerate(self.my_set):
			molecule[0].set_morgan_index(T)
			
'''
txt=open("./data/toxic/corrected_structures/TR000","r").read()
print mol_to_nodes_and_edges(txt)
#a=files_to_nodes_general("./data/mutagenic/mol")
a=files_to_nodes_general("./data/toxic/cancer_nocancer_large/cancer.sdf")
print len(a)
print a[len(a)-1]

b=files_to_nodes_general("./data/toxic/cancer_nocancer_large/nocancer.sdf")
print len(b)
print b[len(b)-1]

'''
'''
b=set_of_molecules()
a=b.smiles_files_to_molecules("./data/cpdb_mutagens/mutagen.tab",1,2)
print b.my_set
#crear una clase molecules set
#creat una clase stats_toolkit
#tratar de arreglar el sobre-entrenamiento

'''
