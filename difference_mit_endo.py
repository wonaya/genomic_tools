import os,sys

### get distribution of differences

mits_val = []
endo_val = []
for a in open("/scratch/02114/wonaya/NCSU_HiSeq/7-21-14_RelicationTiming_Hiseq/merged_1kb/LS_logFC.bedGraph", 'r') :
    mits_val.append(float(a.split("\t")[3].strip("\n")))
for b in open("/scratch/02114/wonaya/NCSU_HiSeq/09-02-14_Maize_endoreplicationTiming_Hiseq/Endo-LS_logFC.bedGraph", 'r') :
    endo_val.append(float(b.split("\t")[3].strip("\n")))
    
diff_val = []
index = 0
for a in mits_val :
    diff_val.append(endo_val[index]-a)
    index += 1

print min(diff_val), max(diff_val)

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

numBins = 50
ax.hist(diff_val,numBins,color='green',alpha=0.8)
#plt.yscale('log',nonposy='clip')
plt.show()
plt.savefig("test.LS.png")

outfile = open("diff.LS.4.bedGraph", 'w')
index = 0
for a in open("/scratch/02114/wonaya/NCSU_HiSeq/7-21-14_RelicationTiming_Hiseq/merged_1kb/LS_logFC.bedGraph", 'r') :
    if float(endo_val[index]-float(a.split("\t")[3].strip("\n"))) >= 4 or float(endo_val[index]-float(a.split("\t")[3].strip("\n"))) <= -4 : 
        outfile.write(a.strip("\n"))
        outfile.write("\t")
        outfile.write(str(endo_val[index]))
        outfile.write("\n")
    index += 1
outfile.close()

### for overlapping region, get gene overlap
os.system("/opt/apps/bedtools/2.19.0/bin/intersectBed -a /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed -b diff.LS.4.bedGraph -wa > diff.LS.4.Zea_mays.AGPv3.23_protein_coding_exon.bed")

### go to mapman annotation or go annotation to get list of terms
genes = []
for a in open("diff.LS.4.Zea_mays.AGPv3.23_protein_coding_exon.bed", 'r') :
    if a.split("\t")[4].split("_")[0] not in genes :
        genes.append(a.split("\t")[4].split("_")[0])


mapman_terms = []
mapman_terms_limited = []
    
for a in open("/work/02114/wonaya/genome/annotation/go_ensembl_zea_mays_annot.txt", 'r') :
    if a.split("\t")[0] in genes :
        mapman_terms.append(a.split("\t")[1].strip("\n"))
        if a.split("\t")[1].strip("\n") not in mapman_terms_limited :
            mapman_terms_limited.append(a.split("\t")[1].strip("\n"))
mapman_terms_limited.sort()
outfile = open("diff.LS.4.Zea_mays.AGPv3.23_protein_coding_exon_go_terms.txt", 'w')
for mapman_term in mapman_terms_limited :
    outfile.write(mapman_term)
    outfile.write("\t")
    outfile.write(str(mapman_terms.count(mapman_term))+"\n")
outfile.close()

