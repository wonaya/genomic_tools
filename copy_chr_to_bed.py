import os,sys

outfile = open("Zea_mays.AGPv3.18_protein_coding_exon_withchr.bed",'w')
for a in open("Zea_mays.AGPv3.18_protein_coding_exon.bed", 'r') :
    outfile.write("chr"+str(a.split("\t")[0])+"\t")
    outfile.write("\t".join(a.split("\t")[1:]))
outfile.close()
