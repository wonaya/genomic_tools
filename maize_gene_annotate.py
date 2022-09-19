import os,sys
import operator
## first make into list non-redundant
nonrd_gene_list = []
for line in open(sys.argv[1], 'r') :
    if line.strip("\n").split("_")[0] not in nonrd_gene_list :
        nonrd_gene_list.append(line.strip("\n").split("_")[0])

term_list = []
term_dict = {}
for go in open("/work/02114/wonaya/genome/annotation/go_ensembl_zea_mays_annot.txt", 'r') :
    if go.split("\t")[0] in nonrd_gene_list and go.split("\t")[2].strip("\n") == "P" :
        term_list.append(go.split("\t")[1])
        term_dict[go.split("\t")[1]] = 0

for term in term_list :
    term_dict[term] += 1

outfile = open(sys.argv[1].strip("_genes.txt")+"_terms.txt", 'w')
for terms in term_dict.keys():
    outfile.write(terms)
    outfile.write("\t")
    outfile.write(str(term_dict[terms]))
    outfile.write("\n")
outfile.close()