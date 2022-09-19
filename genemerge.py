## use python /work/02114/wonaya/scripts/genemerge.py -a gene_goterm.txt -d go_desc.txt -p TAIR10_gene_protein_coding_only.txt -s up.txt 
# export LD_LIBRARY_PATH="/opt/apps/R/2.15.1/lib64/R/lib:$LD_LIBRARY_PATH"
import os,sys
from optparse import OptionParser
import rpy2.robjects as robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
from collections import OrderedDict
 
def extract_gene(infile, name) :
    gene_list = []
    for a in open(infile, 'r') :
        if a.split(",")[0] != '""' :
            if a.split(",")[1].strip('"') not in gene_list :
                gene_list.append(a.split(",")[1].strip('"').split(".")[0])
    outfile = open(name, 'w') 
    for gene in gene_list :
        outfile.write(gene+"\n")
    outfile.close()

def count_gene(assoc, pop, set, desc) :
    import numpy
    desc_dict =get_desc(desc)
    assoc_dict = dict()
    assoc = open(assoc, 'r') 
    assoc_lines = assoc.readlines() 
    for line in assoc_lines :
        for go in line.split("\t")[1].split(";") :
            assoc_dict.setdefault(line.split("\t")[0], []).append(go.strip("\n"))
    print "association dictionary done"
    
    ## dictionary of GOterms P or F
    dict_porf = {}
    for gline in open("gene_association_tair_1.1487_step1", 'r') :
        dict_porf[gline.split("\t")[1]] = gline.split("\t")[3]
    print "porf done"
    pop_count = 0
    pop_list = []
    for gene in open(pop, 'r') :
        pop_count += 1
        pop_list.append(gene.strip("\n"))
    print "pop done"
    if "/" in str(set) :
        setdir = str(set).split("/")[0]
        filename = str(set).split("/")[1]
    else :
        setdir = "."
        filename = str(set)
    exp_evid_code = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"]
    if not os.path.isfile(str(setdir)+"/"+str(filename).strip(".txt")+"_annotation.txt") :
        print "seed annotation file not made"
        sys.exit()
    
    genes_with_exp_evid = []
    rm_redundant_go = {}
    
    
    for aline in open(str(set).strip(".txt")+"_annotation.txt", 'r') :
        rm_redundant_go[aline.split("\t")[1]] = []
        if aline.split("\t")[3].strip("\n") in exp_evid_code :
            if aline.split("\t")[0] not in genes_with_exp_evid :
                genes_with_exp_evid.append(aline.split("\t")[0])
    print "got all GO"
    
    for goterm in rm_redundant_go.keys() :
        for bline in open(str(set).strip(".txt")+"_annotation.txt", 'r') : 
            if bline.split("\t")[1] == goterm and bline.split("\t")[3].strip("\n") in exp_evid_code:
                if bline.split("\t")[0] not in rm_redundant_go[goterm] :
                    rm_redundant_go[goterm].append(bline.split("\t")[0])
    
    print "got all GO term dicts"
    
    set_count = 0
    set_list = []
    for set_gene in open(set, 'r') :
        set_count += 1
        set_list.append(set_gene.strip(" \n"))
    print "got gene lists"
    ## determine percentage of setgenes having at least one experimental evidence
    
    set_go_list = list()
    for assocs in assoc_dict.keys() :
        if assocs in set_list and assoc_dict[assocs] not in set_go_list : 
            set_go_list.extend(assoc_dict[assocs])
    print "got go terms for genes"
    
    set_go_dict = dict()
    for go_ind in set_go_list :
        list_go = []
        for assocs in assoc_dict.keys() :
            if go_ind in assoc_dict[assocs] :
                list_go.append(assocs)
        keys = {}
        for e in list_go  :
            keys[e] = 1
        set_go_dict[go_ind] = list(keys.keys())
        
    print "got set go dict"
    outfile = open(set.split(".txt")[0]+"_hg2_test.out", 'w')
    outfile.write("GO term\tDescription\tP or F\tOverlap Count\tOverlap genelist\tOverlap with Exp.Evid\tGO count\tPopulation\tSet Count\tSet with Exp.Evid\tP-value\tAdj. P-value\n")
    for go_term in set_go_dict.keys() :
        go_count = 0
        go_over_gene_list = []
        for gene_name in set_list :
            if gene_name in set_go_dict[go_term] :
                go_count += 1
                go_over_gene_list.append(gene_name)
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
            robjects.r('adjpval[adjpval>1] <- 1')
            adjpval = robjects.globalenv['adjpval']
            outfile.write(go_term)
            outfile.write("\t")
            outfile.write(desc_dict[go_term])
            outfile.write("\t")
            outfile.write(dict_porf[go_term])
            outfile.write("\t")
            outfile.write(str(r))
            outfile.write("\t")
            outfile.write(",".join(go_over_gene_list))
            outfile.write("\t")
            outfile.write(str(len(rm_redundant_go[go_term])))
            outfile.write("\t")
            outfile.write(str(len(set_go_dict[go_term])))
            outfile.write("\t")
            outfile.write(str(n))
            outfile.write("\t")
            outfile.write(str(k))
            outfile.write("\t")
            outfile.write(str(len(genes_with_exp_evid)))
            outfile.write("\t")
            outfile.write(str(float(str(pval).split(" ")[1].strip("\n"))))
            outfile.write("\t")
            outfile.write(str(float(str(adjpval).split(" ")[1].strip("\n"))))
            outfile.write("\n")
    outfile.close()
    
def get_desc(desc) :
    desc_go = dict()
    desc = open(desc, 'r')
    desc_lines = desc.readlines()
    for desc_line in desc_lines :
        desc_go[desc_line.split("\t")[0]] = desc_line.split("\t")[1].strip("\n")
    return desc_go

def seed_set(assoc, desc, set) :
    outfile = open(str(set).strip(".txt")+"_annotation.txt", 'w')
    assoc_dict = dict()
    assoc = open(assoc, 'r') 
    assoc_lines = assoc.readlines() 
    for line in assoc_lines :
        for go in line.split("\t")[1].split(";") :
            if go.strip("\n") not in assoc_dict.setdefault(line.split("\t")[0], []) :
                assoc_dict.setdefault(line.split("\t")[0], []).append(go.strip("\n"))
    go_dict = get_desc(desc)
    index = 0
    for a in open(set, 'r')  :
        print a.strip("\n"), index
        index += 1
        if a.strip(" \n").upper() in assoc_dict.keys():
            for go in assoc_dict[a.strip(" \n").upper()] :
                for c in open("gene_association_tair_1.1487_step1", 'r') :
                    if c.split("\t")[4].split("|")[0] == a.strip(" \n").upper() and c.split("\t")[1] == go :
                        evid_code = c.split("\t")[2]
                outfile.write(a.strip(" \n").upper())
                outfile.write("\t")
                outfile.write(go)
                outfile.write("\t")
                outfile.write(go_dict[go])
                outfile.write("\t")
                outfile.write(evid_code)
                outfile.write("\n")
    outfile.close()
        
def main():
    parser = OptionParser()
    parser.add_option("-a", dest="assoc", help="association file #80100011\sGO:0009987;GO:0048468;GO:0008219;GO:0012501\n")
    parser.add_option("-d", dest="desc", default="go_desc.txt", help="description file #GO:0000012\ssingle strand break repair\n") 
    parser.add_option("-p", dest="pop", default="TAIR10_gene_protein_coding_only.txt", help="population gene set #80100006\n")
    parser.add_option("-s", dest="set", help="gene set of interest #80100006\n")
    parser.add_option("-i", dest="infile", help="inputfile gene list in csv only valid for 03-24-16")
    (options, args) = parser.parse_args()
    extract_gene(options.infile, options.set)
    seed_set(options.assoc, options.desc, options.set) 
    count_gene(options.assoc, options.pop, options.set, options.desc)
    
if __name__ == '__main__':
    main()