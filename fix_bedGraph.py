import os,sys
import random 
outfile = open(sys.argv[2], 'w') 
for a in open(sys.argv[1], 'r') :
    outfile.write("\t".join(a.split("\t")[:3]))
    outfile.write("\t")
    outfile.write("8C_peaks_"+a.split("\t")[3].strip("\n"))
    outfile.write("\t")
    outfile.write(str(random.randint(1, 10000))+"\t")
    outfile.write(".\t")
    outfile.write(str(random.randint(1, 400))+"\t")
    outfile.write(str(random.randint(1, 400))+"\t")
    outfile.write(str(random.randint(1, 400))+"\n")
outfile.close()