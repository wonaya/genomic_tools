import os,sys

"""
f1 = open("Ath_AGI_LOCUS_TAIR10_Aug2012.txt")
f1lines = f1.readlines()

gene_dict = {}
for f1line in f1lines[1:] :
    if f1line.split("\t")[0] != "'0'" :
        if f1line.split("\t")[2][1:-1].upper() not in gene_dict.keys() :
            if len(f1line.split("\t")[2][1:-1].upper()) > 0 :
                gene_dict[f1line.split("\t")[2][1:-1].upper()] = []

for f1line in f1lines[1:] :
    if f1line.split("\t")[2][1:-1].upper() in gene_dict.keys() :
        gene_dict[f1line.split("\t")[2][1:-1].upper()].append(f1line.split("\t")[0][1:-1])

outfile = open("MapMan_TAIR10_annot.txt", 'w')       
outfile.write("Gene\tAnnot\n")
for gene in gene_dict.keys() :
    for annot in gene_dict[gene] :
        outfile.write(gene)
        outfile.write("\t")
        outfile.write(annot)
        outfile.write("\n")
outfile.close()
"""

f1 = open("ATH_GO_GOSLIM.txt")
f1lines = f1.readlines()

reduced = []
gene_dict = {}
for f1line in f1lines :
    if f1line.split("\t")[7] == 'P' :
        gene_dict[f1line.split("\t")[0]] = []
        reduced.append(f1line)

count = 0
for f1line in reduced :
    if f1line.split("\t")[0] in gene_dict.keys() :
        gene_dict[f1line.split("\t")[0]].append(f1line.split("\t")[5]+"\t"+f1line.split("\t")[3])
    print count, "/", len(reduced)
    count += 1

outfile = open("GOSLIM_20130219_annot.txt", 'w') 
outfile.write("Gene\tAnnot\tRelationship\n")
for gene in gene_dict.keys() :
    for annot in gene_dict[gene] :
        outfile.write(gene)
        outfile.write("\t")
        outfile.write(annot)
        outfile.write("\n")
outfile.close()                
