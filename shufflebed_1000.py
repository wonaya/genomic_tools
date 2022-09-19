### shuffle input bed file 1000 times and calculate overlap

import os,sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy

os.mkdir("tmp")

os.system("/opt/apps/bedtools/2.19.0/bin/intersectBed -wa -a "+sys.argv[1]+" -b /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed > tmp/"+sys.argv[1].split(".bed")[0]+"_overlap_gene.bedGraph")
a = open("tmp/"+sys.argv[1].split(".bed")[0]+"_overlap_gene.bedGraph",'r')
alines = a.readlines()
o = open(sys.argv[1], 'r')
olines = o.readlines()
original_perc = float(len(set(alines)))/float(len(set(olines)))*100
print "no. of peaks original:", len(set(alines)),"/", len(set(olines)),":", original_perc

perm_perc = []
for x in range(0,1000) :
    os.system("/opt/apps/bedtools/2.19.0/bin/shuffleBed -i "+sys.argv[1]+" -g /work/02114/wonaya/genome/maize.genome -excl maize_chrOnly_Ns.bed > tmp/"+str(x)+".bed")
    os.system("/opt/apps/bedtools/2.19.0/bin/intersectBed -wa -a tmp/"+str(x)+".bed -b /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed > tmp/"+str(x).split(".bed")[0]+"_overlap_gene.bedGraph")
    a = open("tmp/"+str(x).split(".bed")[0]+"_overlap_gene.bedGraph", 'r')
    alines = a.readlines()
    perc = float(len(set(alines)))/float(len(set(olines)))*100
    perm_perc.append(perc)

plt.hist(perm_perc, 50, facecolor='green', alpha=0.75)
plt.savefig("perm_1M.png")

os.rmdir("tmp")