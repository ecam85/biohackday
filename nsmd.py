#A set of tools to study Non-sense mediated decay
#Project started at MPS Hack Day 2015 (20/11/2015)

gene_locations = {} #Dictionary {gen_id: (minpos,maxpos)}
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
	

def read_gene_locations(datafile="geneloc.dat"):
	"""
	Reads the gene positions from a data file.
	Format of data file:
	Gene_id\tminpos\tmaxpos
	"""

			


