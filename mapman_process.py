import os,sys

gene_name = []
for a in open(sys.argv[1], 'r') :
    for b in a.split("\r") :
        if b[:2] == "GR" :
            if b[:13] not in gene_name :
                gene_name.append(b[:13])
        else :
            if b.split("_")[0] not in gene_name :
                gene_name.append(b.split("_")[0])

### save mapman term
mapman_term = []
for c in open("Zm_B73_5b_FGS_cds_2012_mapman_summary_gene_list.txt", 'r') :
    if c.split("\t")[1].strip("\n") not in mapman_term :
        mapman_term.append(c.split("\t")[1].strip("\n"))
mapman_term.sort()


mapman_dict = {}
mapman_real_dict = {}
for mapman in mapman_term :
    mapman_dict[mapman] = []
    mapman_real_dict[mapman] = []
    for d in open("Zm_B73_5b_FGS_cds_2012_mapman_summary_gene_list.txt", 'r') :
        if d.split("\t")[1].strip("\n") == mapman :
            mapman_dict[mapman].append(d.split("\t")[0])

for gene in gene_name :
    for mapman in mapman_dict.keys() :
        if gene in mapman_dict[mapman]: 
            mapman_real_dict[mapman].append(gene)
            
        
### save mapman description
mapman_desc = {}
for mapman in mapman_term :
    for e in open("Zm_B73_5b_FGS_cds_2012_mapman_summary_gene_desc.txt", 'r') :
        if e.split("\t")[0] == mapman :
            mapman_desc[mapman] = e.split("\t")[1].strip("\n")

outfile = open("8C_up_MapMan.txt", 'w')
outfile.write("Mapman term\tTerm description\tNo.of set genes in term\tName of set genes in term\n")
for mapman in mapman_term:
    outfile.write(mapman)
    outfile.write("\t")
    outfile.write(mapman_desc[mapman])
    outfile.write("\t")
    outfile.write(str(len(mapman_real_dict[mapman])))
    outfile.write("\t")
    outfile.write(",".join(mapman_real_dict[mapman]))
    outfile.write("\n")

sys.exit()    
"""
### extract top gene from results
outfile3 = open("top98.txt", 'w')
for a in open("top-98-pct.txt", 'r') :
    outfile3.write(a.split(" ")[0])
    outfile3.write("\n")
outfile3.close()
sys.exit()
### extract gene list from GTF
gene_list = []
outfile0 = open("Zea_mays.AGPv2.14_genome.txt", 'w')
for c in open("ccombine_result_filtered.txt", 'r') :
    outfile0.write(c.split("\t")[0])
    outfile0.write("\n")
outfile0.close()

### term \t description \t list of genes
"""
term_list = []
term_gene = {}
for a in open("/work/02114/wonaya/genome/annotation/Zm_B73_5b_FGS_cds_2012.txt" ,'r') :
    if a.split("\t")[0] != "BINCODE" :
        if a.split("\t")[0].strip("'") not in term_list :
            term_list.append(a.split("\t")[0].strip("'"))
            term_gene[a.split("\t")[0].strip("'")] = []

term_desc = {}
for term in term_list :
    for a in open("/work/02114/wonaya/genome/annotation/Zm_B73_5b_FGS_cds_2012.txt" ,'r') :
        if a.split("\t")[0].strip("'") == term :
            term_desc[term] = a.split("\t")[1].strip("'")
            term_gene[term].append(a.split("\t")[2].split("_")[0].upper().strip("'"))
"""
outfile = open("Zm_B73_5b_FGS_cds_2012_mapman_summary_gene_desc.txt", 'w')          
for term in term_list :
    outfile.write(term)
    outfile.write("\t")
    outfile.write(term_desc[term])
    #outfile.write("\t")
    #outfile.write(",".join(term_gene[term]))
    outfile.write("\n")
outfile.close()
"""
outfile2 = open("Zm_B73_5b_FGS_cds_2012_mapman_summary_gene_list.txt", 'w')
for term in term_list :
    for gene in term_gene[term] :
        if len(gene) > 0 :
            outfile2.write(gene)
            outfile2.write("\t")
            outfile2.write(term)
            outfile2.write("\n")
outfile2.close()

