import os,sys

chr = int(sys.argv[1])
outfile = open("CpG_Mo17_all3_bt202.sam_chr"+str(chr)+".methylKit", 'w')
for a in open("CpG_Mo17_all3_bt202.sam_chr678910.methylKit", 'r') :
    if int(a.split("\t")[1]) == chr :
        outfile.write(a)
outfile.close()