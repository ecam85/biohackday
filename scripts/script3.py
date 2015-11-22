#Folows the procedure of script1 (apply all filters and save list) to wt?_chr2L.bam and pcm?_chr2L.bam

#Does not plot.

#Runs from the scripts directory, but output directories must exist.

import os 

os.chdir("..")
print os.getcwd()

import nsmd

pre = "2L/"
list_pcm = ["pcm"+str(i)+"_chr2L.bam" for i in range(1,7)]
list_wt = ["wt"+str(i)+"_chr2L.bam" for i in range(1,7)]
list_out = [pre+"candidates"+str(i)+".txt" for i in range(1,7)]

fccand = pre+"common_candidates.txt"

reg = "2L"

ccand = None

for fpcm,fwt,fout in zip(list_pcm,list_wt,list_out):

	#Minimum reads for PCM: 100
	crit_min=nsmd.crit_min_reads(fpcm,reg,gene=None,minreads=100)

	#Normalized reads in pcm bigger than in wt.
	crit_frac = nsmd.crit_frac_compare(fwt,fpcm,reg)

	#Gene in expressed in pcm with alpha=1
	crit_exp = nsmd.crit_expressed(fpcm,reg)

	#Preliminary candidates satisfying the 3 criteria.
	candidates = [g for g in nsmd.gene_loc if g in crit_min and g in crit_frac and g in crit_exp]

	#Test for center of mass is more expensive. Run it only for candidates.
	crit_cmass = nsmd.crit_cmass(fwt,fpcm,reg,gene=None,only=candidates,alpha=0.1)

	#Save results
	f = open(nsmd.full_path(fout,"results"),"w")

	for g in crit_cmass:
		f.write(g+"\n")
	f.close()

	print fwt+"+"+fpcm+" complete."

