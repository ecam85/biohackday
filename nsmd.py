#A set of tools to study Non-sense mediated decay
#Project started at MPS Hack Day 2015 (20/11/2015)

gene_loc = {} #Dictionary {gen_id: (minpos,maxpos)}
path = {"data":"..data/"} #Paths.

def set_path_data(path,path_type):
	"""
	Sets the path to the data directory.
	"""
	path[path_type]=path

def full_path(name,path_type):
	"""
	Returns the full path to the file "name".
	"""
	return path[path_type]+name
	

def read_gene_loc(datafile="geneloc.dat"):
	"""
	Reads the gene positions from a data file.
	Format of data file:
	Gene_id\tminpos\tmaxpos
	"""
	
	f = open(full_path(datafile,"data"),"r")

	for line in f:
		s = line.rsplit()
		gene_loc[s[0]] = (s[1],s[2])

	f.close()
			
#Gene minloc and maxloc wrappers
def minloc(gene):
	"""
	Returns the min location for the gene.	
	Must be in the gene_loc dictionary.
	"""
	return gene_loc[gene][0]

def maxloc(gene):
	"""
	Returns the max location for the gene.	
	Must be in the gene_loc dictionary.
	"""
	return gene_loc[gene][1]
