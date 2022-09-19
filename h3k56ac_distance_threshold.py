### percentage of ES/MS/LS region in x bp distance of H3K56ac 
import os,sys
import multiprocessing
import readline, glob
import subprocess
from multiprocessing import Process, Manager
"""
### script.py K56_bedgraph Rep_timing_bedGraph

distance_set = 10000 #eg. 10000 for 10kb, 50000 for 50kb

jobs = []

def run_all(chr,return_dict1, return_dict2) :
    total_peaks = 0
    peaks_within_kb = 0
    coord_list = []
    for a in open(sys.argv[1], 'r') :
        if a.split("\t")[0] == str(chr) :
            coord_list.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
    region_list = []
    for b in open(sys.argv[2], 'r') :
        if b.split("\t")[0] == str(chr) : 
            region_list.append(range(int(b.split("\t")[1]),int(b.split("\t")[2])+1))
    count_list = []
    index = 0
    for coord in coord_list :
        count = 0
        for region in region_list :
            ## peaks overlapping
            if len(set(coord)&set(region)) > 0 :
                #print len(set(coord)&set(region)), len(region), len(coord)
                count += 1
            else :
                if coord[-1] < region[0] :
                    if region[0]-coord[-1] < distance_set :
                        count += 1
                elif region[-1] < coord[0] :
                    if coord[0]-region[-1] < distance_set :
                        count += 1    
        count_list.append(count)
        index += 1
    peaks_within_10kb += len(count_list)-count_list.count(0)
    total_peaks += len(count_list)
    return_dict1[chr] = total_peaks
    return_dict2[chr] = peaks_within_kb
    
manager = Manager()
return_dict1 = manager.dict()
return_dict2 = manager.dict()

jobs = []
for chr in range(1,11) :
    print chr
    s1 = multiprocessing.Process(target=run_all, args=(chr,return_dict1, return_dict2))
    jobs.append(s1)
    s1.start()
               
[x.join() for x in jobs]

print sys.argv[1], sys.argv[2], return_dict1.keys(), "total", sum(return_dict1.values()), "peaks within ", distance, "bp", sum(return_dict2.values())

"""

def get_gene(chr,return_dict1,return_dict2,return_dict3,return_dict4) :
    distance_set = 1000 

    ### make file to distinguish 5' or 3' direction
    ### script.py K56_bedgraph 
    gene_direction_dict = {}
    gene_coord_dict = {}
    for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gff3", 'r') :
        if a.split("\t")[0] == str(chr) and a.split("\t")[2] == "gene" :
            gene_direction_dict[a.split("\t")[8].split(";")[0].split(":")[1]] = a.split("\t")[6]
            gene_coord_dict[a.split("\t")[8].split(";")[0].split(":")[1]] = range(int(a.split("\t")[3]),int(a.split("\t")[4])+1)
    
    ### list of K56 peaks
    coord_list = []
    for a in open(sys.argv[1], 'r') :
        if a.split("\t")[0] == str(chr) :
            coord_list.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
    ### how many peaks are within 10kb of gene
    large_list_of_results = []
    print sys.argv[1], "total H3K56ac peak in chr", chr, ":", len(coord_list)
    for coord in coord_list :
        #print chr, coord[0], coord[-1]
        list_of_results = []
        list_of_genes = []
        for gene in gene_coord_dict.keys(): 
            if len(list(set(coord)&set(gene_coord_dict[gene]))) > 0 :
                list_of_results.append("overlap")
                list_of_genes.append(gene)
            elif int(gene_coord_dict[gene][0]) > int(coord[-1]) and int(gene_coord_dict[gene][0])-int(coord[-1]) < distance_set  :
                if gene_direction_dict[gene] == "+" : 
                    list_of_results.append("5")
                    list_of_genes.append(gene)
                    #print chr, gene, gene_coord_dict[gene][0], gene_coord_dict[gene][-1]
                    #sys.exit()
                elif gene_direction_dict[gene] == "-" : 
                    list_of_results.append("3")
                    list_of_genes.append(gene)
            elif int(coord[0]) > int(gene_coord_dict[gene][-1]) and int(coord[0])-int(gene_coord_dict[gene][-1]) < distance_set :
                #print "gene", gene, gene_coord_dict[gene][0], gene_coord_dict[gene][-1]
                #print "peak", coord[0], coord[-1]
                #print gene_direction_dict[gene]
                if gene_direction_dict[gene] == "+" : 
                    list_of_results.append("3")
                    list_of_genes.append(gene)
                elif gene_direction_dict[gene] == "-" : 
                    list_of_results.append("5")
                    list_of_genes.append(gene)
            else :
                list_of_genes.append([])
        if len(list(set(list_of_results))) == 0 :
            large_list_of_results.append(["none"])
        else :
            large_list_of_results.append(list(set(list_of_results)))
        new_list_of_genes = []
        for gene in list_of_genes :
            if len(gene) != 0 :
                new_list_of_genes.append(gene)
        
        #if len(list(set(list_of_results))) != len(new_list_of_genes) :
        #    print chr, coord[0], coord[-1], list(set(list_of_results)), set(new_list_of_genes)
        #    sys.exit()
    print len(large_list_of_results)
    sys.exit()
    count_overlap = 0
    count_5 = 0
    count_3 = 0
    count_none = 0
    for data in large_list_of_results :
        count_overlap += data.count('overlap')
        count_5 += data.count('5')
        count_3 += data.count('3')
        count_none += data.count('none')

    print sys.argv[1], "chr"+str(chr), "overlap", count_overlap, "5'", count_5, "3'", count_3, "none", count_none
    
    return_dict1[chr] = count_overlap
    return_dict2[chr] = count_5
    return_dict3[chr] = count_3
    return_dict4[chr] = count_none
    

def get_genes(chr):
    gene_direction_dict = {}
    gene_coord_dict = {}
    for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gff3", 'r') :
        if a.split("\t")[0] == str(chr) and a.split("\t")[2] == "gene" :
            gene_direction_dict[a.split("\t")[8].split(";")[0].split(":")[1]] = a.split("\t")[6]
            print gene_direction_dict
            gene_coord_dict[a.split("\t")[8].split(";")[0].split(":")[1]] = range(int(a.split("\t")[3]),int(a.split("\t")[4])+1)
            print gene_coord_dict
            sys.exit()
    print len(gene_direction_dict)
    sys.exit()
    ### list of K56 peaks
    coord_list = []
    for a in open(sys.argv[1], 'r') :
        if a.split("\t")[0] == str(chr) :
            coord_list.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
    ### 
    large_list_of_results = []
    print "total gene in chr", chr, ":", len(gene_coord_dict)
    for gene in gene_coord_dict.keys(): 
        print gene
        list_of_results = []
        for coord in coord_list :
            if len(set(coord)&set(gene_coord_dict[gene])) > 0 :
                list_of_results.append("o")
            else :
                list_of_results.append("n")
                """
                if int(gene_coord_dict[gene][0]) > int(coord[-1]) and int(gene_coord_dict[gene][0])-int(coord[-1]) < distance_set :
                    if gene_direction_dict[gene] == "+" : 
                        list_of_results.append("5")
                    elif gene_direction_dict[gene] == "-" : 
                        list_of_results.append("3")
                elif int(coord[0]) > int(gene_coord_dict[gene][-1]) and int(coord[0])-int(gene_coord_dict[gene][-1]) < distance_set :
                    if gene_direction_dict[gene] == "+" : 
                        list_of_results.append("3")
                    elif gene_direction_dict[gene] == "-" : 
                        list_of_results.append("5")
                """
        print list_of_results
        sys.exit()
        #if len(list_of_results) > 1 : 
            #print gene, gene_coord_dict[gene][0], gene_coord_dict[gene][-1], gene_direction_dict[gene], coord[0], coord[-1], list_of_results
            
    
    #count_o = list_of_results.count("o")
    #count_5 = list_of_results.count("5")
    #count_3 = list_of_results.count("3")
    
    #return_dict1[chr] = count_o
    #return_dict2[chr] = count_5
    #return_dict2[chr] = count_3

#get_genes(chr,return_dict1,return_dict2,return_dict3)
#get_genes(chr)
#sys.exit()

    
     
"""
distance_set = 10000 #eg. 10000 for 10kb, 50000 for 50kb

jobs = []

def run_all(chr,return_dict1, return_dict2) :
    total_peaks = 0
    peaks_within_kb = 0
    coord_list = []
    for a in open(sys.argv[1], 'r') :
        if a.split("\t")[0] == str(chr) :
            coord_list.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
    region_list = []
    for b in open(sys.argv[2], 'r') :
        if b.split("\t")[0] == str(chr) : 
            region_list.append(range(int(b.split("\t")[1]),int(b.split("\t")[2])+1))
    count_list = []
    index = 0
    for coord in coord_list :
        count = 0
        for region in region_list :
            ## peaks overlapping
            if len(set(coord)&set(region)) > 0 :
                #print len(set(coord)&set(region)), len(region), len(coord)
                count += 1
            else :
                if coord[-1] < region[0] :
                    if region[0]-coord[-1] < distance_set :
                        count += 1
                elif region[-1] < coord[0] :
                    if coord[0]-region[-1] < distance_set :
                        count += 1    
        count_list.append(count)
        index += 1
    peaks_within_10kb += len(count_list)-count_list.count(0)
    total_peaks += len(count_list)
    return_dict1[chr] = total_peaks
    return_dict2[chr] = peaks_within_kb
    
manager = Manager()
return_dict1 = manager.dict()
return_dict2 = manager.dict()

jobs = []
for chr in range(1,11) :
    print chr
    s1 = multiprocessing.Process(target=run_all, args=(chr,return_dict1, return_dict2))
    jobs.append(s1)
    s1.start()
               
[x.join() for x in jobs]

print sys.argv[1], sys.argv[2], return_dict1.keys(), "total", sum(return_dict1.values()), "peaks within ", distance, "bp", sum(return_dict2.values())
"""

manager = Manager()
return_dict1 = manager.dict()
return_dict2 = manager.dict()
return_dict3 = manager.dict()
return_dict4 = manager.dict()

jobs = []
for chr in range(10,11) :
    s1 = multiprocessing.Process(target=get_gene, args=(chr,return_dict1, return_dict2, return_dict3, return_dict4))
    jobs.append(s1)
    s1.start()
               
[x.join() for x in jobs]

#print sys.argv[1], "total K56 overlap gene", sum(return_dict1.values()), "total K56 at 5'", sum(return_dict2.values()), "total K56 at 3'", sum(return_dict3.values())

### percentage of peaks not overlapping with 
