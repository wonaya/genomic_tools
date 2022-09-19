import os,sys,csv,collections
import numpy
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

def main(count_file):
    base, ext = os.path.splitext(count_file)
    
    print "---->Starting IP1"
    newfile1 = open(sys.argv[1].split(".csv")[0]+"_ip1.txt", 'w')
    for line in open(sys.argv[1], 'r') :
        newfile1.write("\t".join(line.split(",")[:3]))
        newfile1.write("\n")
    newfile1.close()
    counts1 = read_count_file(sys.argv[1].split(".csv")[0]+"_ip1.txt")
    data1, groups1, sizes1, conditions1, genes1 = edger_matrices(counts1)
    probs1 = run_edger(sys.argv[1].split(".csv")[0]+"_ip1.txt")
    outfile1 = "%s_ip1_edgeR_results.txt" % (base)
    write_outfile(outfile1, conditions1, counts1, probs1)
    infile1 = "%s_ip1_edgeR_results.txt" % (base)
    
    filterfile1 = "%s_ip1_edgeR_results_overExpress_all.txt" % (base)
    
    write_filter(infile1, filterfile1)
    
    print "---->Starting IP2"
    newfile2 = open(sys.argv[1].split(".csv")[0]+"_ip2.txt", 'w')
    for line in open(sys.argv[1], 'r') :
        newfile2.write("\t".join(line.split(",")[:2]))
        newfile2.write("\t")
        newfile2.write(line.split(",")[3])
        newfile2.write("\n")
    newfile2.close()
    counts2 = read_count_file(sys.argv[1].split(".csv")[0]+"_ip2.txt")
    data2, groups2, sizes2, conditions2, genes2 = edger_matrices(counts2)
    probs2 = run_edger(sys.argv[1].split(".csv")[0]+"_ip2.txt")
    outfile2 = "%s_ip2_edgeR_results.txt" % (base)
    write_outfile(outfile2, conditions2, counts2, probs2)
    infile2 = "%s_ip2_edgeR_results.txt" % (base)
    
    filterfile2 = "%s_ip2_edgeR_results_overExpress_all.txt" % (base)
    write_filter(infile2, filterfile2)
    os.system("rm -Rf test_edgeR.txt")
    
    overlapfile = "%s_ip1_ip2_overlap_all.txt" % (base)
    txtfile = "%s_annotation_all.txt" % (base)
    overlap(filterfile1, filterfile2, overlapfile,txtfile)
    
def overlap(filterfile1,filterfile2, overlapfile,txtfile):
    outfile = open(overlapfile, 'w')
    txtfile = open(txtfile, 'w')
    outfile.write("Gene\tControl Count\tIP1 count\tIP2 count\tIP1 logFC\tIP2 logFC\tIP1 adj pval\tIP2 adj pval\tGO terms:Description\n")
    list1 = []
    for a in open(filterfile1, 'r') :
        if a.split("\t")[0] != "Gene" :
            list1.append(a.split("\t"))
    list2 = {}
    gene_list = []
    for b in open(filterfile2, 'r') :
        if b.split("\t")[0] != "Gene" :
            list2[b.split("\t")[0]] = b.split("\t")[1:-1]
    for gene1 in list1 :
        if gene1[0] in list2.keys() :
            overlap_list = []
            overlap_list.append(gene1[0])
            overlap_list.append(gene1[1])
            overlap_list.append(gene1[2])
            overlap_list.append(list2[gene1[0]][1])
            overlap_list.append(gene1[4])
            overlap_list.append(list2[gene1[0]][3])
            overlap_list.append(gene1[6])
            overlap_list.append(list2[gene1[0]][5])
            overlap_list.append(gene1[7])
            if str(gene1[7].strip("\n")) != "NA" :
                for x in gene1[7].split(";") :
                    txtfile.write(x.split("=")[1])
                    
            outfile.write("\t".join(overlap_list))
    outfile.close()
    txtfile.close()
    
def write_filter(infile, filterfile):
    ## add annotation 
    ## open annotation file
    term = []
    large_term = {}
    count = 0
    for d in open("annotation/go.obo", 'r') :
        if d[:6] == "[Term]" :
            if count != 0 :
                term2 = []
                for ter in term : 
                    if ter.split(":")[0] == "id" or ter.split(":")[0] == "name" :
                        term2.append(ter.strip("\n"))
                go = term2[0].split(": ")[1]
                desc = term2[1].split(": ")[1]
                large_term[go] = desc  
            term = []
            count += 1
        else : 
            if len(d.split(":")) > 1 :
                term.append(d)
    go_dict = {}
    gene_list = []
    for b in open("annotation/go_ensembl_zea_mays_annot.txt" ,'r') :
        if b.split("\t")[0] != "Gene" and b.split("\t")[2].strip("\n") == "P" : 
            if b.split("\t")[0] not in gene_list : 
                gene_list.append(b.split("\t")[0])
                go_dict[b.split("\t")[0]] = []
    for c in open("annotation/go_ensembl_zea_mays_annot.txt" ,'r') :
        if c.split("\t")[0] != "Gene" and c.split("\t")[2].strip("\n") == "P" : 
            go_dict[c.split("\t")[0]].append(c.split("\t")[1])
    posLogFC = []
    for a in open(infile, 'r') :
        if a.split("\t")[0] == "Gene" :
            labels = []
            for label in a.split("\t") :
                labels.append(label.strip("\r\n"))
            labels.append("GO terms=Description")
            posLogFC.append(labels)
        else:
            #if float(a.split("\t")[4]) > 0 and float(a.split("\t")[6].strip("\n")) <= 0.05 : 
            if float(a.split("\t")[4]) > 0 :
                labels = []
                for label in a.split("\t") :
                    labels.append(label.strip("\r\n"))
                if a.split("\t")[0] in go_dict.keys() :
                    labels_list = []
                    for goterms in  go_dict[a.split("\t")[0]] :
                        if goterms in large_term.keys() :
                            labels_list.append(goterms+"="+large_term[goterms])
                    labels.append(";".join(labels_list))
                else : 
                    labels.append("NA")
                    
                posLogFC.append(labels)
    outfile = open(filterfile, 'w')
    for posLog in posLogFC :
        pos = "\t".join(posLog)
        outfile.write(pos)
        outfile.write("\n")
    outfile.close()    
            
def write_outfile(outfile, conditions, work_counts, probs):
    out_handle = open(outfile, "w")
    writer = csv.writer(out_handle, delimiter='\t')
    writer.writerow(["Gene"] +
            ["%s count" % c for c in conditions] + ["logConc"] + ["logFC"] + ["p-value"] + ["adj.p-value"])
    out_info = []
    for i, gene in enumerate(probs[0]):
        counts = [int(work_counts[c][gene]) for c in conditions]
        out_info.append((probs[1][i], probs[2][i], probs[3][i], probs[4][i], [gene] + counts))
        
    out_info.sort()
    [writer.writerow(start + [prob1] +[prob2]+ [prob3]+ [prob4]) for prob1, prob2, prob3, prob4, start in out_info]

#def run_edger(data, groups, sizes, genes):
def run_edger(in_file):
    """Call edgeR in R and organize the resulting differential expressed genes.
    Changed to estimateCommonDisp from de4DGE (No longer available since bioconductor v.6
    """
    robjects.r('''
        library("edgeR")
    ''')
    # find the version we are running -- check for edgeR exactTest function
    try:
        robjects.r["exactTest"]
        is_13_plus = True
    except LookupError:
        is_13_plus = False
    #params = {'group' : numpy.array(groups), 'lib.size' : sizes}
    #dgelist = robjects.r.DGEList(data, **params)
    os.system("rm -Rf test_edgeR.txt")
    os.system("cp "+str(in_file)+" test_edgeR.txt")
    if is_13_plus:
        # perform Poisson adjustment and assignment as recommended in the manual
        #robjects.globalenv['dP'] = dgelist
        robjects.r('''
            x <- read.delim("test_edgeR.txt", row.names="Gene")
            group <- factor(c(1,2))
            y <- DGEList(counts=x,group=group)
            y <- calcNormFactors(y)
            y <- estimateCommonDisp(y)
            y <- estimateTagwiseDisp(y)
            et <- exactTest(y)
            dP <- topTags(et, n=Inf)
            dP <- data.frame(dP$table)
        ''')
        dgelist = robjects.globalenv['dP']
        indexes = dgelist.rownames
        logconc = list(dgelist[0])
        logfc   = list(dgelist[1])
        pval    = list(dgelist[2])
        adjpval = list(dgelist[3])
    assert len(indexes) == len(pval)
    logconc_w_index = zip(indexes, logconc)
    logconc_w_index.sort()
    logfc_w_index = zip(indexes, logfc)
    logfc_w_index.sort()
    pvals_w_index = zip(indexes, pval)
    pvals_w_index.sort()
    adjpvals_w_index = zip(indexes, adjpval)
    adjpvals_w_index.sort()
    assert len(pvals_w_index) == len(indexes)
    return [k for k,r in logconc_w_index], [r for k,r in logconc_w_index], [s for l,s in logfc_w_index], [q for j,q in pvals_w_index],[p for i,p in adjpvals_w_index]
    #return [r for k,r in logconc_w_index], [s for l,s in logfc_w_index], [q for j,q in pvals_w_index], [p for i,p in adjpvals_w_index]
    
def get_conditions_and_genes(work_counts): 
    conditions = work_counts.keys()
    conditions.sort()
    all_genes = []
    for c in conditions:
        all_genes.extend(work_counts[c].keys())
    all_genes = list(set(all_genes))
    all_genes.sort()
    sizes = []
    for condition in conditions :
        sum = 0
        for gene in all_genes :
            sum += work_counts[condition][gene]
        sizes.append(sum)   
    #sizes = [work_counts[c]["Total"] for c in conditions]
    #print sizes
    #all_genes.remove("Total")
    return conditions, all_genes, sizes
    
def edger_matrices(work_counts):
    """Retrieve matrices for input into edgeR differential expression analysis"""
    conditions, all_genes, sizes = get_conditions_and_genes(work_counts)
    assert len(sizes) == 2
    groups = [1, 2]
    data = []
    final_genes = []
    for g in all_genes:
        cur_row = [int(work_counts[c][g]) for c in conditions]
        if sum(cur_row) > 0:
            data.append(cur_row)
            final_genes.append(g)
    return (numpy.array(data), numpy.array(groups), numpy.array(sizes), conditions, final_genes)

def read_count_file(in_file):
    """Read count information from a simple CSV file into a dictionary"""
    import sys
    counts = collections.defaultdict(dict)
    in_handle = open(in_file, 'r')
    reader = csv.reader(in_handle, delimiter='\t')
    header = reader.next()
    conditions = header[1:3]
    
    for parts in reader:
        region_name = parts[0]
        region_counts = [float(x) for x in parts[1:3]]
        for ci, condition in enumerate(conditions):
            counts[condition][region_name] = region_counts[ci]
    return dict(counts)

if __name__ == "__main__":
    main(sys.argv[1])
