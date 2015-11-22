#Creates a list of candidates satisfying the following critera:
#-Minimum reads in gene, pcm: 100
#-Normalized reads in gene bigger in pcm than in wt.
#-Gene expressed in pcm with alpha=1.
#-Center of mass displaced at least 10% of gene length.
#
#Generates the figures g.jpg and norm_g.jpg, where g is a gene from candidates.

#Data files:
#pcm1_chr2L.bam
#wt1_chr2L.bam

#To run in the same directory as nsmd.py

import nsmd

#Files, region
fpcm = "pcm1_chr2L.bam"
fwt = "wt1_chr2L.bam"
reg = "2L"

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
f = open(nsmd.full_path("candidates1.txt","results"),"w")

for g in crit_cmass:
    f.write(g+"\n")
f.close()

#Save plots:
fignames = [g+".jpg" for g in crit_cmass]
nfignames = ["norm_"+n for n in fignames]

nsmd.plot_pileup(fwt,fpcm,reg,gene=crit_cmass,show=False,filename=fignames)

nsmd.plot_pileup(fwt,fpcm,reg,gene=crit_cmass,show=False,filename=nfignames,norm=True)
