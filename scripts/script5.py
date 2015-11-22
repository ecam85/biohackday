#Order the candidates wrt to number of coincidences in candidates?.txt

import os

os.chdir("..")

import nsmd

pre = "2L/"

list_cand = [pre+"candidates"+str(i)+".txt" for i in range(1,7)]

fout = pre+"candidates_ranking.txt"

c = {}

for fcand in list_cand:
	f = open(nsmd.full_path(fcand,"results"),"r")
	candidates = f.read().splitlines()
	f.close()

	for cc in candidates:
		if cc not in c:
			c[cc] = 1
		else:
			c[cc]+= 1

#Sort by value.
sc = sorted(c,key=c.get)
sc.reverse()

#Save resutls:
f = open(nsmd.full_path(fout,"results"),"w")

for cc in sc:
	f.write(cc+"\t"+str(c[cc])+"\n")

f.close()	


