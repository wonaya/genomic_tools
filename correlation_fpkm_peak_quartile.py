### get top 10 quartile of FPKM 

import os,sys
from scipy import stats
import subprocess
from multiprocessing import Process, Manager
import numpy as np 

list_binary = []
list_cont = []
    
for chr in range(1,11):
    #load in gene and fpkm value as dict
    fpkm_list = []
    for a in open("/scratch/02114/wonaya/NCSU_HiSeq/05-22-13_Maize_Root_RNAseq_HiSeq/RNA-seq_test/merged_cuffout/genes.fpkm_tracking", 'r') :
        if a.split("\t")[6].split(":")[0] == str(chr) and a.split("\t")[0] != "tracking_id" :
            fpkm_list.append(float(a.split("\t")[9]))
    a = np.array(fpkm_list)
    p90 = np.percentile(a, 90)
    
    fpkm_gene_dict = {}
    for a in open("/scratch/02114/wonaya/NCSU_HiSeq/05-22-13_Maize_Root_RNAseq_HiSeq/RNA-seq_test/merged_cuffout/genes.fpkm_tracking", 'r') :
        if a.split("\t")[6].split(":")[0] == str(chr) and a.split("\t")[0] != "tracking_id" and float(a.split("\t")[9]) >= p90 :
            fpkm_gene_dict[a.split("\t")[0]] = a.split("\t")[9]
    
    #load in gene that overlaps with peaks 0 = no overlap, 1 = overlap
    gene_list = []
    gene_coord_dict = {}
    for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gff3", 'r') :
        if a.split("\t")[0] == str(chr) and a.split("\t")[2] == "gene" :
            gene_list.append(a.split("\t")[8].split(";")[0].split(":")[1])
            gene_coord_dict[a.split("\t")[8].split(";")[0].split(":")[1]] = range(int(a.split("\t")[3]),int(a.split("\t")[4])+1)
    
    coord_list = []
    for a in open(sys.argv[1], 'r') :
        if a.split("\t")[0] == str(chr) :
            coord_list.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
    
    list_of_overlap_genes = []
    print "total H3K56ac peak in chr", chr, ":", len(coord_list)
    for coord in coord_list :
        for gene in gene_coord_dict.keys(): 
            if len(list(set(coord)&set(gene_coord_dict[gene]))) > 0 :
                list_of_overlap_genes.append(gene)
    for gene in list(set(gene_list)) :
        if gene in list_of_overlap_genes and gene in fpkm_gene_dict.keys() :
            #print gene, fpkm_gene_dict[gene]
            list_binary.append(1)
            list_cont.append(float(fpkm_gene_dict[gene]))
        elif gene not in list_of_overlap_genes and gene in fpkm_gene_dict.keys() :
            list_binary.append(0)
            list_cont.append(float(fpkm_gene_dict[gene]))
    print "chr"+str(chr), list_binary.count(1), list_binary.count(0)
    
print stats.pointbiserialr(list_binary, list_cont)
print stats.pearsonr(list_binary, list_cont)
print np.corrcoef(list_binary, list_cont)
    

#save into two list binary, continuous
#merge into same array, then do correlation


manager = Manager()
return_dict1 = manager.dict()
return_dict2 = manager.dict()
