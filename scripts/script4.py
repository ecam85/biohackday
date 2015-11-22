#Takes candidadtes lists from script4 and save plots.

import os

os.chdir("..")

import nsmd

pre = "2L/"

list_pcm = ["pcm"+str(i)+"_chr2L.bam" for i in range(1,7)]
list_wt = ["wt"+str(i)+"_chr2L.bam" for i in range(1,7)]
list_cand = [pre+"candidates"+str(i)+".txt" for i in range(1,7)]

list_suffixes = [str(i) for i in range(1,7)]


reg = "2L"

for fcand,suffix,fpcm,fwt in zip(list_cand,list_suffixes,list_pcm,list_wt):
	f = open(nsmd.full_path(fcand,"results"),"r")

	candidates = f.read().splitlines()

	f.close()
	
	if len(candidates) == 0:
		continue	
	
	#Save plots:
	fignames = [pre+g+"_"+suffix+".jpg" for g in candidates]
	nfignames = [pre+g+"_"+suffix+"_norm.jpg" for g in candidates]

	nsmd.plot_pileup(fwt,fpcm,reg,gene=candidates,show=False,filename=fignames)

	nsmd.plot_pileup(fwt,fpcm,reg,gene=candidates,show=False,filename=nfignames,norm=True)
