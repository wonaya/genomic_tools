import os,sys

## gets rid of 4th column
outfile = open(sys.argv[1].split(".")[0]+"_forGMD.bedGraph", 'w')
for a in open(sys.argv[1], 'r') :
    outfile.write(a.split("\t")[0])
    outfile.write("\t")
    outfile.write(a.split("\t")[1])
    outfile.write("\t")
    outfile.write(a.split("\t")[2])
    outfile.write("\t")
    outfile.write(a.split("\t")[4])
    outfile.write("\n")
outfile.close()

#os.system("gzip "+sys.argv[1].split(".")[0]+"_forGAT.bedGraph")