#Generates a plot of fractions of reads per gene for wt1_chr2L.bam and pcm1_chr2L.bam
#Plot name: fracs_wt1_pcm1_2L.jpg
#Must run in same directory as nsmd

import nsmd

#Files, region
fwt = "wt1_chr2L.bam"
fpcm = "pcm1_chr2L.bam"
reg = "2L"
fplot = "fracs_wt1_pcm1_2L.jpg"

#Plot
nsmd.plot_read_comp(fwt,fpcm,reg,show=False,filename=fplot)


