import os,sys


outfile = open(sys.argv[1].split(".wig")[0]+".circos", 'w')
for a in open(sys.argv[1], 'r') :
    outfile.write(a.split("\t")[0])
    outfile.write("\t")
    outfile.write(str((int(a.split("\t")[2])-int(a.split("\t")[1])+1)/2+int(a.split("\t")[1])))
    outfile.write("\t")
    outfile.write(str((int(a.split("\t")[2])-int(a.split("\t")[1])+1)/2+int(a.split("\t")[1])))
    outfile.write("\t")
    outfile.write(a.split("\t")[3])
outfile.close()
