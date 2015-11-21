f = open("got_range.dat","r")
fo = open("geneloc.dat","w")

total = 0
written = 0
ct = 0

for line in f:
	total = total + 1
	s = line.rsplit()
	if len(s) == 3:
		names = s[0].rsplit(",")
		for n in names:
			fo.write(n+"\t"+s[1]+"\t"+s[2]+"\n")
			written = written + 1
	elif len(s)>3:
		#If there is a comma in the first part
		pre = s[0].rsplit(",")
		names = []
		for p in pre:
			names.append(p+"-"+("-".join(s[1:-2])))
		for n in names:
			fo.write(n+"\t"+s[-2]+"\t"+s[-1]+"\n")
			written = written + 1

	else:
		fo.write("NoName"+str(ct)+"\t"+s[-2]+"\t"+s[-1]+"\n")
		ct = ct +1
		written = written + 1

print "Written lines: "+str(written)
print "Total lines: "+str(total)

f.close()
fo.close()


		
