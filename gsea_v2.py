import os,sys

specie = "maize"
## make gene: goterm1, goterm2, goterm3 dictionary
if specie == "arabidopsis" :
    if not os.path.isfile("gene_goterm_C_expevid.txt") :
        gene_dict = {}
        for a in open("gene_association_tair_1.1487_step1", 'r') :
            gene_dict[a.split("\t")[4].split("|")[0].strip("\n")] = []
        exp_code = ["IDA", "IEP", "IMP", "IPI", "IGI"]
        for b in open("gene_association_tair_1.1487_step1", 'r') :
            #if b.split("\t")[3] == "P" and b.split("\t") not in gene_dict[b.split("\t")[4].split("|")[0].strip("\n")]:
            if b.split("\t")[3] == "C" and b.split("\t")[2] in exp_code and b.split("\t") not in gene_dict[b.split("\t")[4].split("|")[0].strip("\n")] :
            #if b.split("\t")[3] == "P" or b.split("\t")[3] == "F" and b.split("\t") not in gene_dict[b.split("\t")[4].split("|")[0].strip("\n")] :
            #if b.split("\t")[3] == "P" or b.split("\t")[3] == "F" and b.split("\t")[2] in exp_code and b.split("\t") not in gene_dict[b.split("\t")[4].split("|")[0].strip("\n")] :
            #if b.split("\t")[3] == "P" and b.split("\t")[2] in exp_code and b.split("\t") not in gene_dict[b.split("\t")[4].split("|")[0].strip("\n")] :
                gene_dict[b.split("\t")[4].split("|")[0].strip("\n")].append(b.split("\t")[1])
        outfile1 = open("gene_goterm_C_expevid.txt", 'w')
        for gene in gene_dict.keys() :
            if len(gene_dict[gene]) > 0 :
                outfile1.write(gene+"\t"+(";").join(gene_dict[gene])+"\n")
        outfile1.close()
    
    # make go : description file
    if not os.path.isfile("go_desc.txt") :
        c1 = open("gene_ontology_ext.obo", 'r')
        c1lines = c1.readlines()
        go_desc = []
        for c1line in c1lines :
            if c1line[:2] == "id" :
                #go_desc[c1line[3:].strip(" ").strip("\n")] = c1lines[c1lines.index(c1line)+1][6:].strip("\n")
                go_desc.append(c1line[3:].strip(" ").strip("\n")+"\t"+c1lines[c1lines.index(c1line)+1][6:].strip("\n"))
        outfile2 = open("go_desc.txt", 'w')
        for go in go_desc :
            outfile2.write(go+"\n")
        outfile2.close()
    # save all protein coding gene into a file
    outfile3 = open("TAIR10_gene_protein_coding_only.txt" , 'w') 
    for d in open("TAIR10_gene_type", 'r') :
        if d.split("\t")[1].strip("\n") == "protein_coding" :
            outfile3.write(d.split("\t")[0]+"\n")
    outfile3.close()
    
    outfile4 = open("go_summary.txt", 'w')
    all_go = []
    for e in open("gene_goterm_P.txt", 'r') :
        if len(e.split("\t")) > 1 :
            for f in e.split("\t")[1].split(";") :
                if f.strip("\n") not in all_go :
                    all_go.append(f.strip("\n"))
    for go in all_go :
        count = 0
        gene_list = []
        for g in open("go_desc.txt", 'r') :
            if g.split("\t")[0] == go :
                desc = g.split("\t")[1].strip("\n")
        for h in open("gene_goterm.txt", 'r') :
            if len(h.split("\t")) == 2 :
                for i in h.split("\t")[1].split(";") :
                    if i.strip("\n") == go :
                        count += 1
                        if h.split("\t")[0] not in gene_list :
                            gene_list.append(h.split("\t")[0])
        if len(go) == 10 :
            outfile4.write(go)
            outfile4.write("\t")
            outfile4.write(desc)
            outfile4.write("\t")
            outfile4.write(str(count))
            outfile4.write("\t")
            outfile4.write(",".join(gene_list))
            outfile4.write("\n")       
    outfile4.close()
    
    # run gene_merge.pl
    set_file = sys.argv[1] ### list of genes
    hg_file = set_file.split(".txt")[0]+"_hg"
    hg_outfile = set_file.split(".txt")[0]+"_hg_out.txt"
    sort_hg_outfile = set_file.split(".txt")[0]+"_hg_sort.txt"
    os.system("/work/02114/wonaya/scripts/genemerge.pl -a ../gene_goterm_PF_allevid.txt -d ../go_desc.txt -p ../TAIR10_gene_protein_coding_only.txt -s "+set_file+" -o "+hg_file)


elif specie == "maize" :
    if not os.path.isfile("maize_gene_mapman.txt") :
        gene_dict = {}
        for a in open("/work/02114/wonaya/genome/annotation/MapMan_B73_annot.txt", 'r') :
            if a.split("\t")[0] != "Gene" :
                gene_dict[a.split("\t")[0]] = []
        for a in open("/work/02114/wonaya/genome/annotation/MapMan_B73_annot.txt", 'r') :
            if a.split("\t")[0] != "Gene" :
                gene_dict[a.split("\t")[0]].append(a.split("\t")[1].strip("\n"))
        outfile1 = open("maize_gene_mapman.txt", 'w')
        for gene in gene_dict.keys() :
            if len(gene_dict[gene]) > 0 :
                outfile1.write(gene+"\t"+(";").join(gene_dict[gene])+"\n")
        outfile1.close()
        
    if not os.path.isfile("mapman_desc.txt") :
        outfile = open("mapman_desc.txt", 'w') 
        for a in open("/work/02114/wonaya/genome/annotation/Zm_B73_5b_FGS_cds_2012_mapman_summary.txt", 'r') :
            outfile.write("\t".join(a.split("\t")[:2]))
            outfile.write("\n")
        outfile.close()
    
    if not os.path.isfile("maize_protein_coding_only.txt") :
        outfile = open("maize_protein_coding_only.txt", 'w') 
        gene = []
        for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon.bed", 'r') :
            gene.append(a.split("\t")[3].split("_")[0])
        set_gene = list(set(gene))
        for gene in set_gene :
            outfile.write(gene)
            outfile.write("\n")
        outfile.close()
    # run gene_merge.pl
    set_file = sys.argv[1] ### list of genes
    hg_file = set_file.split(".")[0]+"_hg"
    os.system("perl /work/02114/wonaya/scripts/genemerge_v2.pl -a maize_gene_mapman.txt -d mapman_desc.txt -p maize_protein_coding_only.txt -s "+set_file+" -o "+hg_file)

hg_file = set_file.split(".")[0]+"_hg.txt"
outfile = open(set_file.split(".")[0]+"_hg_rmNA.txt", 'w')
for a in open(hg_file, 'r') :
    if a.split("\t")[7] != "NA" :
        outfile.write(a)
outfile.close()

# sort
pval_dict = {}
pval_list = []
for f in open(set_file.split(".")[0]+"_hg_rmNA.txt", 'r') :
    pval_dict[float(f.split("\t")[7])] = []
    if float(f.split("\t")[7]) not in pval_list :
        pval_list.append(float(f.split("\t")[7]))

for g in open(set_file.split(".")[0]+"_hg_rmNA.txt", 'r') :
    tmp_list = []
    for gi in g.split("\t")[:7] :
        tmp_list.append(gi)
    tmp_list.append(g.split("\t")[8].strip("\n"))
    pval_dict[float(g.split("\t")[7].strip("\n"))].append(tmp_list)
pval_list.sort()

outfile5 = open(set_file.split(".")[0]+"_hg_sort.txt", 'w') 
outfile5.write("term\tdescription\ttotal genes in term\ttotal protein coding genes\ttotal seed genes in term\ttotal seed genes\tp-val\tbonferroni corrected p-val\n")
for pval in pval_list :
    for val in pval_dict[pval] :
        outfile5.write(val[0])
        outfile5.write("\t")
        outfile5.write(val[-1])
        outfile5.write("\t")
        outfile5.write("\t".join(val[2:6]))
        outfile5.write("\t")
        outfile5.write(str(val[-2])+"\t")
        outfile5.write("\t")
        outfile5.write(str(pval)+"\n")
outfile5.close()
#os.system("rm "+hg_outfile)
