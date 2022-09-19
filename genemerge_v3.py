import os,sys
import numpy
from scipy import stats
import rpy2.robjects as robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr

### hg analysis
try :
    specie = sys.argv[1]
except IndexError:
    print "please provide a specie, arabidopsis or maize"
    sys.exit()

### step 1. make association file gene_goterm
if specie == "maize" :
    desc_go = dict()
    desc = open("/work/02114/wonaya/genome/annotation/go_desc_maize.txt", 'r')
    desc_lines = desc.readlines()
    for desc_line in desc_lines :
        desc_go[desc_line.split("\t")[0]] = desc_line.split("\t")[1].strip("\n")
    
    pop_list = []
    for gene in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed", 'r') :
        pop_list.append(gene.split("\t")[3].split("_")[0].strip("\n"))
    pop_list_2 = list(set(pop_list))
    pop_count = len(pop_list_2)
    pop_list = pop_list_2
    
    ### seed set 
    set_list = []
    for set_gene in open(sys.argv[2], 'r') :
        set_list.append(set_gene.strip(" \n"))
    set_list_2 = list(set(set_list))
    set_list = set_list_2
    set_count = len(set_list)
    
    ### gene: go dict
    assoc_dict = dict()
    assoc = open("/work/02114/wonaya/genome/annotation/maize_gene_mapman.txt", 'r') 
    assoc_lines = assoc.readlines() 
    for line in assoc_lines :
        for go in line.split("\t")[1].split(";") :
            assoc_dict.setdefault(line.split("\t")[0], []).append(go.strip("\n"))
    
    set_go_list = list()
    for gene in assoc_dict.keys() :
        if gene in set_list and assoc_dict[gene] not in set_go_list : 
            set_go_list.extend(assoc_dict[gene])
    set_go_list2 = list(set(set_go_list))
    set_go_list = set_go_list2
    
    ### go: gene dict
    set_go_dict = dict()
    for go_term in set_go_list :
        for gene in assoc_dict.keys() :
            if go_term in assoc_dict[gene] :
                if gene not in set_go_dict.setdefault(go_term, []) :
                    set_go_dict.setdefault(go_term, []).append(gene)
    
    outfile = open(sys.argv[2].split(".")[0]+"_hg.out", 'w')
    outfile.write("GO term\tDescription\tOverlap Count\tGO count\tPopulation\tSet Count\tP-value\tAdj. P-value\n")
    for go_term in set_go_dict.keys() :
        go_count = 0
        for gene_name in set_list :
            if gene_name in set_go_dict[go_term] :
                go_count += 1
        n = pop_count
        p = go_count/float(pop_count)
        k = len(set_list)
        r = go_count
        if go_count > 1 :
            #print go_term, go_count-1, len(set_go_dict[go_term]), n-(len(set_go_dict[go_term])), k
            q = 1-p
            np = int(n*p+0.5)
            nq = int(n*q+0.5)
            x = n
            sum = 0
            while x > n-k :
                sum += numpy.log(int(x))
                x -= 1
            y = 2
            sum_2 = 0
            while y <= k :
                sum_2 += numpy.log(int(y))
                y += 1
            i = n*(1-p)
            sum_3 = 0
            while i > n*(1-p)-(k-np) :
                sum_3 += numpy.log(int(i))
                i -= 1
            j = 2
            sum_4 = 0
            while j <= (k-np) :
                sum_4 += numpy.log(int(j))
                j += 1
            
            robjects.r('count <- '+str(r-1))
            robjects.r('r <- '+str(r))
            robjects.r('gocount <- '+str(len(set_go_dict[go_term])))
            robjects.r('pop <- '+str(int(n)-len(set_go_dict[go_term])))
            robjects.r('set <- '+str(k))
            robjects.r('pval <- phyper(count,gocount,pop,set,lower.tail=FALSE)')
            pval = robjects.globalenv['pval']
            robjects.r('adjpval <- phyper(count,gocount,pop,set,lower.tail=FALSE)*r')
            adjpval = robjects.globalenv['adjpval']
            outfile.write(go_term)
            outfile.write("\t")
            outfile.write(desc_go[go_term])
            outfile.write("\t")
            outfile.write(str(r))
            outfile.write("\t")
            outfile.write(str(len(set_go_dict[go_term])))
            outfile.write("\t")
            outfile.write(str(n))
            outfile.write("\t")
            outfile.write(str(k))
            outfile.write("\t")
            outfile.write(str(float(str(pval).split(" ")[1].strip("\n"))))
            outfile.write("\t")
            outfile.write(str(float(str(adjpval).split(" ")[1].strip("\n"))))
            outfile.write("\n")
    outfile.close()
    
