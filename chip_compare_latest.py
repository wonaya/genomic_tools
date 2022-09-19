#!/usr/bin/env python
  
import os,sys
import numpy as np
import matplotlib.mlab as mlab
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

### column dependent ### get rid of this

def permutation(static_set, test_set, coordinate_N, permutation_cycle, bedtools_dir, unique_out, fai) :
    try :
        os.mkdir("tmp")
    except OSError :
        os.chdir(".")
    
    os.system(bedtools_dir+"/intersectBed -wa -a "+test_set+" -b "+static_set+" > tmp/"+test_set.split(".")[0]+"_overlap_"+static_set.split(".")[0]+".bedGraph")
    a = open("tmp/"+test_set.split(".")[0]+"_overlap_"+static_set.split(".")[0]+".bedGraph",'r')
    alines = a.readlines()
    o = open(test_set, 'r')
    olines = o.readlines()
    original_perc = float(len(set(alines)))/float(len(set(olines)))*100
    print "no. of peaks original:", len(set(alines)),"/", len(set(olines)),":", original_perc
    
    print "\nstart time:", datetime.datetime.now()
    perm_perc = []
    for x in range(0,1000) :
        os.system(bedtools_dir+"/shuffleBed -i "+test_set+" -g /work/02114/wonaya/genome/maize.genome -excl "+coordinate_N+" > tmp/"+str(x)+".bed")
        os.system(bedtools_dir+"/intersectBed -wa -a tmp/"+str(x)+".bed -b "+static_set+" > tmp/"+str(x).split(".")[0]+"_overlap_"+static_set.split(".")[0]+".bedGraph")
        a = open("tmp/"+str(x).split(".bed")[0]+"_overlap_"+static_set.split(".")[0]+".bedGraph", 'r')
        alines = a.readlines()
        perc = float(len(set(alines)))/float(len(set(olines)))*100
        perm_perc.append(perc)
    print "end time:", datetime.datetime.now() # 5 min for 1000 peaks 1000 times, 
    
    count = 0
    for perm in perm_perc :
        if perm >= original_perc : 
            count += 1
    print "p-value", float(count)/float(1000)*100
    plt.hist(perm_perc, 50, facecolor='green', alpha=0.75)
    plt.xlabel('percentage overlap (%)')
    plt.ylabel('frequency')
    plt.plot([original_perc, original_perc], [0,100], 'k-', lw=2, color='r')
    plt.savefig("perm_1K.png")
    
    if unique_out != None :
        whole_static_file = open(static_set, 'r') 
        whole_static_lines = whole_static_file.readline()
        static_column_count = whole_static_lines.split("\t")
        whole_test_file = open(test_set, 'r')
        whole_test_lines = whole_test_file.readline()
        test_column_count = whole_test_lines.split("\t")
        print static_column_count, test_column_count
    
        chr_list = []
        static_dict_1 = {}
        static_dict_2 = {}
        for a in open(fai, 'r') :
            chr_list.append(a.split("\t")[0])
            static_dict_1[a.split("\t")[0]] = []
            static_dict_2[a.split("\t")[0]] = []
        ### generate unique set for static set
        list_1 = []
        for a in open(static_set, 'r') :
            if a.split("\t")[0] in static_dict_1.keys() :
                static_dict_1[a.split("\t")[0]].append(a.split("\t")[1]+"-"+a.split("\t")[2])
        list_2 = []
        for a in open(test_set, 'r') :
            if a.split("\t")[0] in static_dict_2.keys() :
                static_dict_2[a.split("\t")[0]].append(a.split("\t")[1]+"-"+a.split("\t")[2])
        outfile = open(static_set.split(".")[0]+"_unique.bedGraph", 'w')
        for chr in chr_list :
            unique_peaks_test = []
            uniq_count = 0
            for a in static_dict_1[chr] : 
                for b in static_dict_2[chr] : 
                    print a, b
                    if len(set(range(int(a.split("-")[0]),int(a.split("-")[1])+1))&set(range(int(b.split("-")[0]),int(b.split("-")[1])+1))) > 0 :
                        uniq_count += 1
                if uniq_count == 0 :
                    unique_peaks_test.append(a)
            print chr, unique_peaks_test
        
    
def jaccard(test_set, static_set) :
    jaccard_list = []
    for chr in range(1,11) :
        list_a = []
        list_b = []
        for a in open(test_set, 'r') :
            if a.split("\t")[0] == str(chr) :
                list_a.extend(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
        for b in open(static_set, 'r') :
            if b.split("\t")[0] == str(chr) :
                list_b.extend(range(int(b.split("\t")[1]),int(b.split("\t")[2])+1))
        intersection_cardinality = len(set.intersection(*[set(list_a), set(list_b)]))
        union_cardinality = len(set.union(*[set(list_a), set(list_b)]))
        jaccard_list.append(intersection_cardinality/float(union_cardinality))
    print np.mean(jaccard_list)

def point_biserial(chip, rnaseq, min_percentage_overlap,gff,offset,bedtools_dir,macs_score) :
    ### how many overlap over 50%
    ### gene list, gene : fpkm
    genes = [] 
    gene_fpkm_dict = {}
    for a in open(rnaseq, 'r') :
        if a.split("\t")[0] != "tracking_id" :
            genes.append(a.split("\t")[0])
            gene_fpkm_dict[a.split("\t")[0]] = float(a.split("\t")[9])
    genes_bed = open("tmp/gene.bed", 'w')
    gene_coord_dict = {}
    for a in open(gff, 'r') :
        if a[0] != "#" and a.split("\t")[2] == "gene" :
            gene = a.split("\t")[8].split(";")[0].split(":")[1]
            gene_coord_dict[gene] = a.split("\t")[3]+"-"+a.split("\t")[4]
            genes_bed.write(str(a.split("\t")[0]+"\t"+a.split("\t")[3]+"\t"+a.split("\t")[4]+"\t"+gene+"\n"))
    genes_bed.close()
    filter_file = open("tmp/chip.bed", 'w')
    for a in open(chip, 'r') :
        if float(a.split("\t")[6]) > float(macs_score) :
            filter_file.write(a)
    filter_file.close()
    os.system(bedtools_dir+"/intersectBed -wao -a tmp/chip.bed -b tmp/gene.bed > tmp/"+chip.split(".")[0]+"_overlap_"+rnaseq.split(".")[0]+".bedGraph")
    ### 
    overlap_gene = []
    for a in open("tmp/"+chip.split(".")[0]+"_overlap_"+rnaseq.split(".")[0]+".bedGraph",'r') :
        if a.split("\t")[13].strip("\n") != "." :
            overlap = a.split("\t")[13].strip("\n")
            peak_distance = int(a.split("\t")[2])-int(a.split("\t")[1])
            if float(overlap)/float(peak_distance)*100 >= float(min_percentage_overlap) :
                overlap_gene.append(a.split("\t")[12])
    
    chip_list = []
    rna_list = []
    gene_list = []
    for gene in genes :
        ### does it have overlap?
        if gene in overlap_gene:
            chip_list.append(1)
            rna_list.append(math.log(gene_fpkm_dict[gene]+1,2))
            gene_list.append(gene)
        else :
            chip_list.append(0)
            rna_list.append(math.log(gene_fpkm_dict[gene]+1,2))
            gene_list.append(gene)
    rvalue = stats.pointbiserialr(chip_list, rna_list)[0] 
    print "r-value", rvalue
    
def ttest(chip, rnaseq, min_percentage_overlap,gff,offset,bedtools_dir,percentile,macs_score_percentile) :
    ### how many overlap over 50%
    ### gene list, gene : fpkm
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
    percentile_cutoff = np.percentile(fpkm_list, percentile)
    
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
    print "MACs score top", macs_score_percentile, "% =", np.percentile(chip_val_list, float(macs_score_percentile))
    filter_file.close()
    os.system(bedtools_dir+"/intersectBed -wao -a tmp/chip.bed -b tmp/gene.bed > tmp/"+chip.split(".")[0]+"_overlap_"+rnaseq.split(".")[0]+".bedGraph")
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
    print "p-value for test :", p_val, "p-value for random :", random_p_val
    
    bins = np.linspace(-20, 20, 100)

    ax1 = plt.subplot(211)
    plt.hist(t_test_1, bins, alpha=0.5, label='Peak', color='red')
    plt.hist(t_test_0, bins, alpha=0.5, label='No Peak', color='blue')
    plt.ylabel('Log(frequency)')
    plt.legend(loc='upper right')
    ax2 = plt.subplot(212)
    plt.hist(r_test_1, bins, alpha=0.5, label='Peak', color='red')
    plt.hist(r_test_0, bins, alpha=0.5, label='No Peak', color='blue')
    
    plt.legend(loc='upper right')
    plt.xlabel('Log(FPKM)')
    plt.ylabel('Log(frequency)')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    plt.show()
    
def bisttest(chip, bsseq, bedtools_dir) :
    try :
        os.mkdir("tmp")
    except OSError :
        os.chdir(".")
    genes = [] 
    for a in open(bsseq, 'r') :
        print a
    percentile_cutoff = np.percentile(fpkm_list, percentile)
    

def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma

def rnaseq_moving_average(chip, rnaseq, min_percentage_overlap,gff,offset,bedtools_dir) :
    ### how many overlap over 50%
    ### gene list, gene : fpkm
    genes = [] 
    gene_fpkm_dict = {}
    for a in open(rnaseq, 'r') :
        if a.split("\t")[0] != "tracking_id" :
            genes.append(a.split("\t")[0])
            gene_fpkm_dict[a.split("\t")[0]] = float(a.split("\t")[9])
    genes_bed = open("tmp/gene.bed", 'w')
    gene_coord_dict = {}
    for a in open(gff, 'r') :
        if a[0] != "#" and a.split("\t")[2] == "gene" :
            gene = a.split("\t")[8].split(";")[0].split(":")[1]
            gene_coord_dict[gene] = a.split("\t")[3]+"-"+a.split("\t")[4]
            genes_bed.write(str(a.split("\t")[0]+"\t"+a.split("\t")[3]+"\t"+a.split("\t")[4]+"\t"+gene+"\n"))
            #genes_bed.write(str(gene_fpkm_dict[gene])+"\n")
    genes_bed.close()
    os.system(bedtools_dir+"/intersectBed -wa -a tmp/gene.bed -b "+chip+" > tmp/"+chip.split(".")[0]+"_overlap_"+rnaseq.split(".")[0]+".bedGraph")
    overlap_gene_list = []
    for a in open("tmp/"+chip.split(".")[0]+"_overlap_"+rnaseq.split(".")[0]+".bedGraph" ,'r') :
        overlap_gene_list.append(a.split("\t")[3].strip("\n"))
    os.system("sort -nrk4 tmp/gene_fpkm.bed > tmp/gene_fpkm_sorted.bed")
    outfile = open("moving_average.txt", 'w')
    for a in open("tmp/gene_fpkm_sorted.bed", 'r'):
        if float(a.split("\t")[3]) > 0 :
            outfile.write(a.split("\t")[4].strip("\n"))
            outfile.write("\t")
            outfile.write(a.split("\t")[3])
            outfile.write("\t")
            if a.split("\t")[4].strip("\n") in overlap_gene_list :
                outfile.write("1\n")
            else : 
                outfile.write("0\n")
    outfile.close()
    
    gene_label = []
    data = []
    count_list = []
    count = 1
    for a in open("moving_average.txt", 'r') :
        gene_label.append(a.split("\t")[0])
        data.append(float(a.split("\t")[1].strip("\n")))
        count_list.append(count)
        count += 1
    column_labels = gene_label
    row_labels = ['FPKM']
    #data_np = np.array(count_list,data)
    plt.subplots()
    ax1 = plt.subplot(2,1,1)
    ax1.plot(count_list,data)
    ax1.set_ylabel('log(FPKM)')
    ax1.set_xlabel('gene')
    #my_data = np.genfromtxt("moving_average.txt", delimiter="\t")
    #heatmap = ax1.pcolor(my_data,vmin=0, vmax=10000)
    #ax1.set_xticks(np.arange(data_np.shape[1])+0.5, minor=False)
    #ax1.set_yticks(np.arange(data_np.shape[0])+0.5, minor=False)
    #ax1.invert_yaxis()
    #ax1.xaxis.tick_top()
    #ax1.set_xticklabels(column_labels, minor=False)
    #ax1.set_yticklabels(column_labels, minor=False)
    #plt.colorbar(heatmap)
    x = []
    y = []
    for a in open("moving_average.txt",'r') :
        x.append(float(a.split("\t")[1]))
        y.append(float(a.split("\t")[2].strip("\n"))) 
 
    yMA = movingaverage(y,200)
    ax2 = plt.subplot(2,1,2)
    ax2.plot(x[len(x)-len(yMA):],yMA)
    #ax2.set_xlim([0.5, 1] )
    #ax2.invert_yaxis()
    plt.show()
    

parser = OptionParser()
allgroup = OptionGroup(parser, "Required for all function")
allgroup.add_option("--type", dest="type", help="experiment type ChIP/ChIP or ChIP/RNAseq or ChIP/BSseq")
allgroup.add_option("--test_1", dest="static", help="not permutated set")
allgroup.add_option("--test_2", dest="test", help="permutated set/RNAseq/BSseq")
allgroup.add_option("--bedtools_dir", dest="bedtools_dir", help="full path of bedtools binary")
parser.add_option_group(allgroup)
prmgroup = OptionGroup(parser, "Required for permutation")
prmgroup.add_option("--Ncoord", dest="ncoord", help="coordinates of N in bedgraph format /opt/apps/bedtools/2.22.1/bin/ --Ncoord /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/maize_chrOnly_Ns.bed")
prmgroup.add_option("--permutate_n", dest="permutate_n", help="no_of_permutation_cycles_to_run", default=1000)
prmgroup.add_option("--unique_out", action="store_false", dest="unique_out", help="output unique peaks as separate bedGraph file")
prmgroup.add_option("--fai", dest="fai", help="genome index file")
parser.add_option_group(prmgroup)
rnagroup = OptionGroup(parser, "Required for ChIP/RNAseq")
rnagroup.add_option("--rnaseq", dest="rnaseq", help="genes.fpkm_tracking from cufflinks")
rnagroup.add_option("--gff3", dest="gff3", help="gff3 file")
rnagroup.add_option("--minoverlap", dest="minoverlap", help="minimum percentage overlap between ChIP/gene (% of peak)", default=0)
rnagroup.add_option("--offset", dest="offset", help="distance in bp from +/- gene to be tested", default=0)
rnagroup.add_option("--percentile", dest="percentile", help="percentile of FPKM values to be tested", default=0)
rnagroup.add_option("--macs_score", dest="macs_score", help="threshold of macs score percentile", default=0)
parser.add_option_group(rnagroup)
bsgroup = OptionGroup(parser, "Required for ChIP/BSseq")
bsgroup.add_option("--bsseq", dest="bsseq", help="bedgraph from Bisulfite pipeline")
parser.add_option_group(bsgroup)
(options, args) = parser.parse_args()
   
def main():
    if options.type == "ChIP/ChIP" :
        permutation(options.static,options.test,options.ncoord,options.permutate_n,options.bedtools_dir,options.unique_out,options.fai) 
        #point_biserial("2C_peaks.broadPeak","genes.fpkm_tracking", 0, "/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gff3",0,"/opt/apps/bedtools/2.22.1/bin")
    elif options.type == "ChIP/ChIP jac" :
        jaccard(options.static,options.test)
    elif options.type == "ChIP/RNAseq" :
        ttest(options.static,options.rnaseq, options.minoverlap, options.gff3,options.offset,options.bedtools_dir,options.percentile,options.macs_score)
        point_biserial(options.static,options.rnaseq, options.minoverlap, options.gff3,options.offset,options.bedtools_dir,options.macs_score)
        
if __name__ == "__main__":
    main()   
