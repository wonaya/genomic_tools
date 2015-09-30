#!/usr/bin/env python
  
import os,sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
import datetime
from math import *
from scipy import stats
from optparse import OptionParser, OptionGroup
import pandas as pd
from urllib2 import urlopen
import math
import random

#def ttest(chip, rnaseq, gff, macs_score_percentile, min_percentage_overlap) :
def ttest(chip, rnaseq, gff, macs_score_percentile) :
    min_percentage_overlap = 0
    print chip, rnaseq, gff, macs_score_percentile
    try :
        os.mkdir("tmp")
    except OSError :
        os.chdir(".")
    
    genes = [] 
    gene_fpkm_dict = {}
    random_dict = {}
    fpkm_list = []
    for a in open(rnaseq, 'r') :
        if a.split("\t")[0] != "tracking_id" :
            genes.append(a.split("\t")[0])
            gene_fpkm_dict[a.split("\t")[0]] = float(a.split("\t")[9])
            random_dict[a.split("\t")[0]] = random.uniform(0, 10000)
            fpkm_list.append(float(a.split("\t")[9]))
    
    #percentile_cutoff = np.percentile(fpkm_list, percentile)
    
    genes_bed = open("tmp/gene.bed", 'w')
    gene_coord_dict = {}
    for a in open(gff, 'r') :
        if a[0] != "#" and a.split("\t")[2] == "gene" and a.split("\t")[8].split(";")[0].split(":")[1] in gene_fpkm_dict.keys():
            gene = a.split("\t")[8].split(";")[0].split(":")[1]
            gene_coord_dict[gene] = a.split("\t")[3]+"-"+a.split("\t")[4]
            genes_bed.write(str(a.split("\t")[0]+"\t"+a.split("\t")[3]+"\t"+a.split("\t")[4]+"\t"+gene+"\t"))
            genes_bed.write(str(gene_fpkm_dict[gene])+"\n")
    genes_bed.close()
    chip_val_list = []
    for a in open(chip, 'r') :
        chip_val_list.append(float(a.split("\t")[6]))
    filter_file = open("tmp/chip.bed", 'w')
    for a in open(chip, 'r') :
        if float(a.split("\t")[6]) > np.percentile(chip_val_list, float(macs_score_percentile)) :
            filter_file.write(a)
    #print "MACs score top", macs_score_percentile, "% =", np.percentile(chip_val_list, float(macs_score_percentile))
    filter_file.close()
    os.system("intersectBed -wao -a tmp/chip.bed -b tmp/gene.bed > tmp/"+chip.split(".")[0]+"_overlap_"+rnaseq.split(".")[0]+".bedGraph")
    overlap_gene = []
    for a in open("tmp/"+chip.split(".")[0]+"_overlap_"+rnaseq.split(".")[0]+".bedGraph",'r') :
        if a.split("\t")[13].strip("\n") != "." :
            overlap = a.split("\t")[14].strip("\n")
            peak_distance = int(a.split("\t")[2])-int(a.split("\t")[1])
            if float(overlap)/float(peak_distance)*100 >= float(min_percentage_overlap) :
                overlap_gene.append(a.split("\t")[12])
    
    t_test_1 = []
    t_test_0 = []
    r_test_1 = []
    r_test_0 = []
    for gene in genes :
        ### does it have overlap?
        if gene in overlap_gene :
            if float(gene_fpkm_dict[gene]) > 0 :
                t_test_1.append(math.log(gene_fpkm_dict[gene],2))
                r_test_1.append(math.log(random_dict[gene],2))
        else :
            if float(gene_fpkm_dict[gene]) > 0 :
                t_test_0.append(math.log(gene_fpkm_dict[gene],2))
                r_test_0.append(math.log(random_dict[gene],2))
    t_stat, p_val = stats.ttest_ind(t_test_0, t_test_1, equal_var=False)
    r_stat, random_p_val = stats.ttest_ind(r_test_0, r_test_1, equal_var=False)
    #print "p-value for test :", p_val, "p-value for random :", random_p_val
    
    outfile = open(chip.split(".")[0]+"_overlap_"+rnaseq.split(".")[0]+".txt", 'w') 
    outfile.write("MACs score top "+str(macs_score_percentile)+"% = "+str(np.percentile(chip_val_list, float(macs_score_percentile)))+"\n")
    outfile.write("p-value for test : "+str(p_val)+" p-value for random : "+str(random_p_val)+"\n")
    outfile.close()
    bins = np.linspace(-20, 20, 100)

    ax1 = plt.subplot(111)
    plt.hist(t_test_1, bins, alpha=0.5, label='Peak', color='red')
    plt.hist(t_test_0, bins, alpha=0.5, label='No Peak', color='blue')
    plt.xlabel('Log(FPKM)')
    plt.ylabel('Log(frequency)')
    plt.legend(loc='upper right')
    ax1.set_yscale('log')
    ax1.set_yscale('log')
    plt.savefig(chip.split(".")[0]+"_overlap_"+rnaseq.split(".")[0]+"_test.png")
    plt.clf()
    
    ax2 = plt.subplot(111)
    plt.hist(r_test_1, bins, alpha=0.5, label='Peak', color='red')
    plt.hist(r_test_0, bins, alpha=0.5, label='No Peak', color='blue')
    plt.legend(loc='upper right')
    plt.xlabel('Log(FPKM)')
    plt.ylabel('Log(frequency)')
    ax2.set_yscale('log')
    ax2.set_yscale('log')
    plt.savefig(chip.split(".")[0]+"_overlap_"+rnaseq.split(".")[0]+"_random.png")
    os.system("rm -Rf tmp")

parser = OptionParser()
rnagroup = OptionGroup(parser, "Required for ChIP/RNAseq")
rnagroup.add_option("--chip", dest="static", help=".broadPeak output from MACS2 --broad")
rnagroup.add_option("--rnaseq", dest="rnaseq", help="genes.fpkm_tracking from cufflinks")
rnagroup.add_option("--gff3", dest="gff3", help="gff3 file")
rnagroup.add_option("--macs_score", dest="macs_score", help="threshold of macs score percentile", default=0)
#rnagroup.add_option("--minoverlap", dest="minoverlap", help="minimum percentage overlap between ChIP/gene (% of peak)", default=0)
parser.add_option_group(rnagroup)
(options, args) = parser.parse_args()

def main():
    #ttest(options.static,options.rnaseq,options.gff3,options.macs_score,options.minoverlap)
    ttest(options.static,options.rnaseq,options.gff3,options.macs_score)     

if __name__ == "__main__":
    main()   
