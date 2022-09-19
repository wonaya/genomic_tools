import os,sys

outfile = open(sys.argv[2], 'w')
a = open(sys.argv[1],'r')
alines = a.readlines()
for aline in alines[0].split("\r"):
    if aline[:3] == "chr" :
        outfile.write(aline.split("\t")[0][3:])
        outfile.write("\t")
        outfile.write("\t".join(aline.split("\t")[1:4]))
        outfile.write("\n")
outfile.close()
sys.exit()
    

