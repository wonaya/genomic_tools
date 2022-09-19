import os,sys

count = 0
for a in open("/work/02114/wonaya/genome/annotation/Zea_mays.AGPv3.18_protein_coding_exon.bed", 'r') :
    count += int(a.split("\t")[2])-(int(a.split("\t")[1])-1)
print count