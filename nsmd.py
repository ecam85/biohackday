#A set of tools to study Non-sense mediated decay
#Project started at MPS Hack Day 2015 (20/11/2015)

#Authors:
#	Miroslav Batchkarov
#	Nick Ayres
#	Eduard Campillo-Funollet
#	Heather McAslan
#	Amrita Ganpatial

#To-Do:
#	Get gene information from genome file.

import pysam

gene_loc = {} #Dictionary {gen_id: (minpos,maxpos)}
path = {"data":"../data/"} #Paths.

#Runs read_gene_loc to get locations from file.
read_gene_loc()

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
		gene_loc[s[0]] = (int(s[1]),int(s[2]))

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

#Read counting.
def read_count(bamfile,region,gene=None):
	"""
	Count total number of reads in a gene or in the whole region.
	Note that the gene must be in the specified region and in the geneloc dict.
	"""

	samfile = pysam.Samfile(full_path(bamfile,"data")) #Open the Bamfile.	
	if not gene:
		return samfile.count(region)
	else:
		return samfile.count(region,minloc(gene),maxloc(gene))

def gene_read_frac(bamfile,region,gene,cached_count=None):
	"""
	Total reads in a gene over total reads in a region.
	If cached_count is given, it is used as total reads in region.
	"""
	if not cached_count:
		cached_count = read_count(bamfile,region)

	return float(read_count(bamfile,region,gene))/cached_count
 
