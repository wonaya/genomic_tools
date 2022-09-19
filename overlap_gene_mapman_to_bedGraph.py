import os,sys

enriched_terms = []
for a in open("0-5mm_enriched_genes_hg_sort.txt", 'r') :
    if a.split("\t")[0] != "term" :
        if float(a.split("\t")[8].strip("\n")) <= 0.05 :
            enriched_terms.append(a.split("\t")[0])

all_gene_in_term = []
for a in open("maize_gene_mapman.txt", 'r') :
    if a.split("\t")[1].strip("\n") in enriched_terms :
        all_gene_in_term.append(a.split("\t")[0])

gene_in_term_enriched = []    
for a in open("0-5mm_enriched_genes.txt", 'r') :
    if a.strip("\n") in all_gene_in_term :
        gene_in_term_enriched.append(a.strip("\n"))

### get coordinates for these genes

outfile = open("0-5mm_sig_enriched_genes.bed", 'w')
for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gff3", 'r') :
    if a[0] != "#" :
        if a.split("\t")[2] == "gene" :
            if a.split("\t")[8].split(";")[0].split(":")[1] in gene_in_term_enriched : 
                outfile.write(a.split("\t")[0])
                outfile.write("\t")
                outfile.write(a.split("\t")[3])
                outfile.write("\t")
                outfile.write(a.split("\t")[4])
                outfile.write("\t")
                outfile.write(a.split("\t")[8].split(";")[0].split(":")[1])
                outfile.write("\n")
outfile.close()

os.system("sort -k1 -k2 0-5mm_sig_enriched_genes.bed > 0-5mm_sig_enriched_genes.sorted.bed")