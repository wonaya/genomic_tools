#### converts fai to bed 
#### output chr	start	end	chrom_name

import os, sys
outfile = open(sys.argv[1].split(".fai")[0]+".bed",'w')
for line in open(sys.argv[1], 'r') :
    outfile.write(line.split("\t")[0])
    outfile.write("\t")
    outfile.write("0\t")
    outfile.write(line.split("\t")[1])
    outfile.write("\t")
    outfile.write("chrom\n")
outfile.close()
