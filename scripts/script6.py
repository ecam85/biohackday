#Plot CG7294 for all the samples.

list_pcm = ["pcm"+str(i)+"_chr2L.bam" for i in range(1,7)]
list_wt = ["wt"+str(i)+"_chr2L.bam" for i in range(1,7)]
list_out = ["CG7294_"+str(i)+".jpg" for i in range(1,7)]
list_nout = ["CG7294_"+str(i)+"_norm.jpg" for i in range(1,7)]

reg = "2L"
g = "CG7294"

import os

os.chdir("..")

import nsmd

for fpcm,fwt,fout,fnout in zip(list_pcm,list_wt,list_out,list_nout):
	nsmd.plot_pileup(fwt,fpcm,reg,g,False,fout)
	nsmd.plot_pileup(fwt,fpcm,reg,g,False,fnout,True)


	


