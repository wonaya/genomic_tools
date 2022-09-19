#USE python fai windowsize out

import os,sys

outfile = open(sys.argv[3], 'w')
windowsize = int(sys.argv[2])
for a in open(sys.argv[1], 'r') :
    final = float(a.split("\t")[1])/float(windowsize)
    abs= str(final).split(".")[0]
    remainder = float(a.split("\t")[1])-int(abs)*float(windowsize)
    chrom = a.split("\t")[0]
    #print a.split("\t")[0], a.split("\t")[1], final, remainder
    for x in range(0, int(abs)+1) :
        outfile.write(a.split("\t")[0]+"\t")
        outfile.write(str(x*windowsize)+"\t"+str((x+1)*windowsize)+"\n")
    outfile.write(a.split("\t")[0]+"\t"+str(x*windowsize)+"\t"+str((x*windowsize)+remainder)+"\n")
outfile.close()
