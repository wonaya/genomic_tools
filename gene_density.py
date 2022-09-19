import os,sys
### calculate gene density based on gff3 
"""
outfile = open("Zea_mays.AGPv3.23_gene.bed", 'w')
for chr in range(1,11) :
    a = open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gff3", 'r')
    alines = a.readlines()
    for aline in alines[1:] :
        if aline.split("\t")[0] == str(chr) and aline.split("\t")[2] == "gene" :
            outfile.write(aline.split("\t")[0])
            outfile.write("\t")
            outfile.write(aline.split("\t")[3])
            outfile.write("\t")
            outfile.write(aline.split("\t")[4])
            outfile.write("\n")
outfile.close()

max_dict = {}
for a in open("/work/02114/wonaya/genome/maize.genome" ,'r') :
    max_dict[a.split("\t")[0]] = a.split("\t")[1].strip("\n")

import math
def roundup(x):
    return int(math.ceil(x / 10000)+1) * 10000
outfile = open("chr_10kb_maize.bed", 'w')
for chr in range(1,11) :
    for x in range(0,(roundup(int(max_dict[str(chr)]))/10000)+1) :
        outfile.write(str(chr))
        outfile.write("\t")
        outfile.write(str(x*10000+1))
        outfile.write("\t")
        outfile.write(str((x+1)*10000))
        outfile.write("\n")
outfile.close()
"""

### draw histogram of gene density 
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

x = []
for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/chr_10kb_maize_gene_density.bedGraph", 'r') :
    x.append(float(a.split("\t")[3].strip("\n")))
plt.hist(x, 50, facecolor='green', alpha=0.75)
plt.show()

