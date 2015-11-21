#A set of tools to study Non-sense mediated decay
#Project started at MPS Hack Day 2015 (20/11/2015)

#Authors:
#	Miroslav Batchkarov
#	Nick Ayres
#	Eduard Campillo-Funollet
#	Heather McAslan
#	Amrita Ganpatial

#To-Do:
#	Gene location file including region (chromosome)
#	Get gene information from genome file.
#	Global cached file counts.
#	Pysam version control. Now based in Pysam 0.6
#	Criterium wrapper - all crit_* functions share structure.

import pysam #SAM/BAM files handling.
from matplotlib import pyplot as plt #Plotting

gene_loc = {} #Dictionary {gen_id: (minpos,maxpos)}
path = {"data":"./data/","fig":"./fig/"} #Paths.

def set_path(path,path_type):
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
	
#Runs read_gene_loc to get locations from file.
read_gene_loc()

		
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

def get_region_length(bamfile,region):
	"""
	Returns the length of the region
	"""
	sf = pysam.Samfile(full_path(bamfile))
	ret = sf.lengths[sf.gettid(region)]
	sf.close()

	return ret

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

	samfile.close()



def gene_read_frac(bamfile,region,gene=None,cached_count=None):
	"""
	Total reads in a gene over total reads in a region.
	If cached_count is given, it is used as total reads in region.
	If not gene, a mapping for all in gene_loc.
	"""
	if not cached_count:
		cached_count = read_count(bamfile,region)
	
	if gene:
		return float(read_count(bamfile,region,gene))/cached_count

	gene_frac = {}

	for gene in gene_loc:
		gene_frac[gene] = gene_read_frac(bamfile,region,gene,cached_count)

	return gene_frac

def gene_read_total(bamfile,region,gene=None):
	"""
	Total reads in a gene.
	If not gene, a mapping for all in gene_loc.
	"""
	
	if gene:
		return read_count(bamfile,region,gene)

	ret = {}

	for gene in gene_loc:
		ret[gene] = gene_read_total(bamfile,region,gene)

	return ret 


def gene_reads_compare(file1,file2,region,gene=None,ccount1=None,ccount2=None):
	"""
	Compares the proportion of reads in a gene from two different BAM files.

	If not gene, returns a dictionary {gene:(frac1,frac2)} for all genes in geneloc.

	If gene, returns the pair (frac1,frac2) for the gene.

	Cached counts used if passed.
	"""

	if not ccount1:
		ccount1 = read_count(file1,region)

	if not ccount2:
		ccount2 = read_count(file2,region)

	if gene:
		return (gene_read_frac(file1,region,gene,cached_count=ccount1),gene_read_frac(file2,region,gene,cached_count=ccount2))

	#If gene not passed, do it for all genes.
	ret = {}
	for gene in gene_loc:
		ret[gene] = gene_reads_compare(file1,file2,region,gene,ccount1,ccount2) 

	return ret

def plot_read_comp(file1,file2,region):
	"""
	Given two files an a region, plots for each gene the fraction of total reads in that gene"

	Red +  for the fraction in file 1, blue o for the fraction in file 2.

	Note: the gene is plotted in the mid-point of each location.
	"""

	read_frac = gene_reads_compare(file1,file2,region)

	xx = []
	yy1 = []
	yy2 = []

	for gene in gene_loc:
		xx.append((maxloc(gene)+minloc(gene))/2.)
		yy1.append(read_frac[gene][0])
		yy2.append(read_frac[gene][1])

	plt.figure()
	plt.plot(xx,yy1,"+",color="red")
	plt.plot(xx,yy2,"o",color="blue")
	plt.show(block=False)
	
def crit_frac_compare(file1,file2,region,gene=None,only=None,ccount1=None,ccount2=None):
	"""
	If not gene,
	Returns a list of genes in gene_loc such that

	The fraction (wrt total reads in region) of reads in the gene in file1 is smaller or equal than the fraction in file2.

	Note: we expect that candidate genes in PCM will have in general more reads than in WT, once normalized wrt to total reads.
	If gene, returns true if the criteria is true for that gene.	
	If only, consider only genes in the list "only".
	"""

	if not ccount1:
		ccount1 = read_count(file1,region)

	if not ccount2:
		ccount2 = read_count(file2,region)

	if gene:
		return gene_read_frac(file1,region,gene,cached_count=ccount1) <= gene_read_frac(file2,region,gene,cached_count=ccount2)

	#If gene not passed, do it for all genes.
	ret = [] 
	for gene in gene_loc:
		if not only or gene in only:		
			if crit_frac_compare(file1,file2,region,gene,ccount1=ccount1,ccount2=ccount2):
				ret.append(gene)

	return ret

def crit_min_reads(bamfile,region,gene=None,minreads=0,frac=False,only=None,ccount=None):
	"""
	Returns true if gene satisfies the following criteria.

	The total number of reads is bigger or equal than minreads.

	If frac, fractional number of reads over total in region is used.

	If not gene, a returns a list of all the genes in gene_loc that satisfy the criterium. If only, only genes in only.
	"""

	if not gene and frac and not ccount:
		ccount = read_count(bamfile,region)

	if gene:
		if frac:
			return gene_read_frac(bamfile,region,gene,ccount) >= minreads

		else:
			return gene_read_total(bamfile,region,gene) >= minreads

	ret = []

	for gene in gene_loc:
		if not only or gene in only:
			if crit_min_reads(bamfile,region,gene,minreads,frac,only,ccount):
				ret.append(gene)

	return ret

def is_expressed(bamfile,region,gene,region_length=None,alpha = 1.0,ccount=None):
	"""
	Returns True if gene satisfies:
	
	Reads in gene/(gene length) >= alpha* reads in region/(region length)
	"""

	if not ccount:
		ccount = read_count(bamfile,region)
	
	if not region_length:
		region_length = get_region_length(bamfile,region)

	reads_in_gene = read_count(bamfile,region,gene)
	gene_length = maxloc(gene)-minloc(gene)
	

	return alpha*ccount/region_length <= float(reads_in_gene)/gene_length

def crit_expressed(bamfile,region,gene=None,only=None,region_length=None,alpha=1.0,ccount=None):
	"""
	Returns True if gene is expressed with factor alpha.

	If not gene, for all in gene_loc. If only, only if in only.
	"""

	if not ccount:
		ccount = read_count(bamfile,region)

	if not region_length:
		region_length=get_region_length(bamfile,region)

	if gene:
		return is_expressed(bamfile,region,gene,region_length,alpha,ccount)

	ret = []
	ct = 0
	for gene in gene_loc:
		if not only or gene in only:
			if crit_expressed(bamfile,region,gene,only,region_length,alpha,ccount):
				ret.append(gene)

	return ret

#def get_pileup(bamfile,region,gene=None):
#	"""
#	NOTE: This is not doing what we expected!
#	Returns the counts for each location for all the reads that cover that location. 
#	"""
#	samfile = pysam.Samfile(full_path(bamfile,"data"))
#	if gene:
#		pileup = samfile.pileup(region,minloc(gene),maxloc(gene))
#	else:
#		pileup = samfile.pileup(region)
#
#	ret = []
#
#	for read in pileup:
#		ret.append(read.n)
#
#	samfile.close()
#	return ret

def get_pileup(bamfile,region,gene=None):
	"""
	Returns the counts for each location for all the reads that cover that location. 
	Pileup by hand, not using pileup function.
	"""
	samfile = pysam.Samfile(full_path(bamfile,"data"))
	if gene:
		samiter = samfile.fetch(region,minloc(gene),maxloc(gene))
		ret=[0]*(maxloc(gene)-minloc(gene)+200) #200 for safety reasons.
		shift = minloc(gene)
	else:
		samiter = samfile.fetch(region)
		ret = [0]*(get_region_length(bamfile,region)+200)
		shift = 0

	p = []
	l = []

	for read in samiter:
		p.append(read.pos)
		l.append(read.rlen)	

	for i in range(len(p)):
		for j in range(l[i]):
			ret[p[i]-shift+j] += 1	

	samfile.close()
	return ret

def cmass(bamfile,region,gene,norm=False, ccount=None):
	"""
	Center of mass of the gene normalized pileup wrt total reads in region.

	If norm, center of mass in [0,1] with respect to gene length.
	"""

	if not ccount:
		ccount = read_count(bamfile,region)
	
	pileup = get_pileup(bamfile,region,gene)
	norm_pileup = [float(p)/ccount for p in pileup]

	s = sum(norm_pileup)

	if s==0:
		if norm:
			return .5
		else:
			return (minloc(gene)+maxloc(gene))/2.
		 

	shift = minloc(gene) #pileup starts at 0, but minloc!=0 in general.

	c = sum([loc*norm_pileup[loc-shift] for loc in range(minloc(gene),shift+len(norm_pileup))])/sum(norm_pileup)

	if not norm:
		return c
	else:
		return (c-shift)/(maxloc(gene)-minloc(gene))

def cmass_diff(file1,file2,region,gene,norm=False,ccount1=None,ccount2=None):
	"""
	Return the difference of centers of mass between gene in file1 and gene in file2, c2-c1

	if norm, normalized to 0,1 wrt gene length.
	""" 		

	if not ccount1:
		ccount1 = read_count(file1,region)

	if not ccount2:
		ccount2 = read_count(file2,region)

	c1 = cmass(file1,region,gene,norm,ccount1)
	c2 = cmass(file2,region,gene,norm,ccount2)

	return c2-c1

def crit_cmass(file1,file2,region,gene=None,only=None,alpha=.1,ccount1=None,ccount2=None):
	"""
	Returns the gene satisfying the following criterium:

	True if the normalized center of mass is displaced more than alpha.
	If not gene, for all in gene_loc. If only, only genes in only.
	"""
	
	if not ccount1:
		ccount1 = read_count(file1,region)

	if not ccount2:
		ccount2 = read_count(file2,region)

	if gene:
		return abs(cmass_diff(file1,file2,region,gene,True,ccount1,ccount2)) >= alpha

	ret = []
	for gene in gene_loc:
		if not only or gene in only:
			if crit_cmass(file1,file2,region,gene,only,alpha,ccount1,ccount2):
				ret.append(gene)
	
	return ret


def plot_pileup(file1,file2=None,region="2L",gene=None,show=True,filename="test.jpg"): 
	"""
	Plots the pileup histogram for the given files.
	If file2, for both files (file1 in red, file2 in blue).
	If not gene, for the whole region.
	If gene is a list, a figure for each gene.
	
	If show=False, save the result to test.jpg in the fig path.
	If show=False and gene is a list, filename must be a list.
	"""

	pileup2 = None

	if isinstance(filename,str):
		filename = [filename]

	#Get pileups.
	if not gene:
		pileup1 = get_pileup(file1,region)
		if file2:
			pileup2 = get_pileup(file2,region)
	elif isinstance(gene,str):
		pileup1 = get_pileup(file1,region,gene)
		if file2:
			pileup2=get_pileup(file2,region,gene)
	else:
		pileup1 = [get_pileup(file1,region,g) for g in gene]
		if file2:
			pileup2 = [get_pileup(file2,region,g) for g in gene]

	#Plot

	nfig = []
	if not gene or isinstance(gene,str):
		fig=plt.figure()
		plt.plot(pileup1,color="red")
		if pileup2:
			plt.plot(pileup2,color="blue")
		nfig.append(fig.number)
	else:
		for i in range(len(gene)):
			fig=plt.figure()
			plt.plot(pileup1,color="red")
			if pileup2:
				plt.plot(pileup2,color="blue")
			nfig.append(fig.number)

	#Show or save
	if show:
		plt.show(block=False)
	else:
		for i in range(len(gene)):
			plt.figure(nfig[i])
			fig.savefig(full_path(filename[i],"fig"))


			
		
	
		
