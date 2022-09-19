import os,sys

gene_list = []
for a in open(sys.argv[1],'r') :
    gene_list.append(a.split("\t")[3].strip("\n").split(".")[0])

gene_list = list(set(gene_list))
outfile = open(sys.argv[1].split(".bed")[0]+"_gene_list.txt", 'w')
for gene in gene_list :
    outfile.write(gene)
    outfile.write("\n")
outfile.close()