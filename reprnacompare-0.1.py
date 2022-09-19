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

def boxplot(reptiming, rnaseq, gff,max_fpkm,min_fpkm, min_percentage_overlap, name) :
    try :
        os.mkdir("tmp")
    except OSError :
        os.chdir(".")
    
    ### filter FPKM file
    fpkm_list = []
    for a in open(rnaseq, 'r') :
        if a.split("\t")[0] != "tracking_id" :
                fpkm_list.append(float(a.split("\t")[9]))
    filter_file = open("tmp/fpkm.out", 'w')
    for a in open(rnaseq, 'r') :
        if a.split("\t")[0] != "tracking_id" :
            if float(a.split("\t")[9]) >= np.percentile(fpkm_list, float(min_fpkm)) and float(a.split("\t")[9]) <= np.percentile(fpkm_list, float(max_fpkm)):
                filter_file.write(a)
    filter_file.close()
    
    ### save all FPKM value to a dictionary
    genes = [] 
    gene_fpkm_dict = {}
    random_dict = {}
    for a in open("tmp/fpkm.out", 'r') :
        genes.append(a.split("\t")[0])
        gene_fpkm_dict[a.split("\t")[0]] = float(a.split("\t")[9])
    
    ### gene coordinates
    genes_bed = open("tmp/gene.bed", 'w')
    gene_coord_dict = {}
    for a in open(gff, 'r') :
        if a[0] != "#" and a.split("\t")[2] == "gene" and a.split("\t")[8].split(";")[0].split(":")[1] in gene_fpkm_dict.keys():
            gene = a.split("\t")[8].split(";")[0].split(":")[1]
            gene_coord_dict[gene] = a.split("\t")[3]+"-"+a.split("\t")[4]
            genes_bed.write(str(a.split("\t")[0]+"\t"+a.split("\t")[3]+"\t"+a.split("\t")[4]+"\t"+gene+"\t"))
            genes_bed.write(str(gene_fpkm_dict[gene])+"\n")
    genes_bed.close()
    
    ### get list of different types of segment
    reptiming_dict = {}
    reptiming_val_dict = {}
    reptiming_count_dict = {}
    reptiming_max_dict = {}
    reptiming_min_dict = {}
    reptiming_box_dict = {}
    for a in open(reptiming, 'r') :
        if a[0] != "#" :
            reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
            reptiming_val_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
            reptiming_count_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
            reptiming_max_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
            reptiming_min_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
            reptiming_box_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
    for a in open(reptiming, 'r') :
        if a[0] != "#" :
            chr = a.split("\t")[0]
            coord_start = int(a.split("\t")[3])-1
            coord_end   = int(a.split("\t")[4])-1
            reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]].append(chr+":"+str(coord_start)+"-"+str(coord_end))
    
    #txtfile = open("reptiming_overlap_"+rnaseq.split(".")[0]+".txt", 'w') 
    
    ### write into bed and sort
    print "\nstart time:", datetime.datetime.now()
    for time in reptiming_dict.keys() :
        outfile = open("tmp/"+time+".bed", 'w') 
        for coord in reptiming_dict[time] :
            outfile.write(coord.split(":")[0])
            outfile.write("\t")
            outfile.write(coord.split(":")[1].split("-")[0])
            outfile.write("\t")
            outfile.write(coord.split(":")[1].split("-")[1])
            outfile.write("\t"+str(time)+"\n")
        outfile.close()
        #os.system("sortBed -i tmp/"+time+".bed > tmp/"+time+"_sorted.bed")
    
        os.system("intersectBed -wao -a tmp/"+time+".bed -b tmp/gene.bed > tmp/"+time+"_overlap_"+rnaseq.split(".")[0]+".bedGraph")
        overlap_gene = []
        for a in open("tmp/"+time+"_overlap_"+rnaseq.split(".")[0]+".bedGraph",'r') :
            if a.split("\t")[4] != "." :
                #overlap_gene.append(a.split("\t")[7])
                overlap = a.split("\t")[-1].strip("\n")
                peak_distance = int(a.split("\t")[2])-int(a.split("\t")[1])
                if float(overlap)/float(peak_distance)*100 >= float(min_percentage_overlap) :
                    overlap_gene.append(a.split("\t")[-3])
        overlap = list(set(overlap_gene))
        ### write overlap gene names
        genelist_file = open(time+"_overlap_gene_list.txt", 'w')
        for overlap_gene in overlap :
            genelist_file.write(overlap_gene+"\n")
        genelist_file.close()
        count_0 = 0
        count_above = 0
        for gene in genes :
            ### does it have overlap? only above 0
            if gene in overlap :
                if float(gene_fpkm_dict[gene]) > 0 :
                     #reptiming_val_dict[time].append(gene_fpkm_dict[gene])
                     reptiming_val_dict[time].append(math.log(gene_fpkm_dict[gene],2))
                     count_above += 1
                elif float(gene_fpkm_dict[gene]) == 0 :
                     count_0 += 1
        
        #print time, count_0, count_above
        B= plt.boxplot(reptiming_val_dict[time])
        reptiming_box_dict[time].append(str(time)+"\t"+str([item.get_ydata()[1] for item in B['whiskers']][0])+"\t"+str([item.get_ydata()[0] for item in B['whiskers']][0])+"\t"+str(np.median(reptiming_val_dict[time]))+"\t"+str([item.get_ydata()[0] for item in B['whiskers']][1])+"\t"+str([item.get_ydata()[1] for item in B['whiskers']][1]))
        
        reptiming_max_dict[time] += count_above
        reptiming_min_dict[time] += count_0
    
    max_val = [reptiming_max_dict['ES'],reptiming_max_dict['ESMS'],reptiming_max_dict['MS'],reptiming_max_dict['MSLS'],reptiming_max_dict['LS'],reptiming_max_dict['ESLS'],reptiming_max_dict['ESMSLS']]
    min_val = [reptiming_min_dict['ES'],reptiming_min_dict['ESMS'],reptiming_min_dict['MS'],reptiming_min_dict['MSLS'],reptiming_min_dict['LS'],reptiming_min_dict['ESLS'],reptiming_min_dict['ESMSLS']]
    
    ind = np.arange(len(reptiming_max_dict))   
    width= 0.35
    fig,  ax= plt.subplots()
    rects1 = ax.bar(ind, max_val, width,color='r')
    rects2 = ax.bar(ind+width, min_val, width,color='y')
    ax.set_ylabel('Gene Count')
    ax.set_xticks(ind+width)
    ax.set_xticklabels(('ES','ESMS','MS','MSLS','LS','ESLS','ESMSLS'))
    ax.legend((rects1[0], rects2[0]), ('>0', '=0'))
    
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,'%d' % int(height),ha='center', va='bottom')
    autolabel(rects1)
    autolabel(rects2)
    plt.savefig(name+'_barchart.png')
    plt.close()
    for reptiming in reptiming_val_dict.keys(): 
        reptiming_count_dict[reptiming] += len(reptiming_val_dict[reptiming])
    es_test = np.array(reptiming_val_dict['ES'])
    esms_test = np.array(reptiming_val_dict['ESMS'])
    esls_test = np.array(reptiming_val_dict['ESLS'])
    esmsls_test = np.array(reptiming_val_dict['ESMSLS'])
    ms_test = np.array(reptiming_val_dict['MS'])
    msls_test = np.array(reptiming_val_dict['MSLS'])
    ls_test = np.array(reptiming_val_dict['LS'])
    data = [es_test,esms_test,ms_test,msls_test,ls_test,esls_test,esmsls_test]
    fig = plt.figure()
    ax1= fig.add_subplot(1,1,1)
    ax1.boxplot(data,showfliers=False)
    #ax1.boxplot(data)
    ax1.set_ylabel('log(FPKM)')
    ax1.set_xticklabels(['ES', 'ESMS','MS', 'MSLS', 'LS','ESLS', 'ESMSLS'])
    ax1.get_xaxis().tick_bottom()
    plt.savefig(name+'_boxplot.png')
    boxfile = open(name+"_boxplot_value.txt", 'w')
    for time in ['ES', 'ESMS','MS', 'MSLS', 'LS','ESLS', 'ESMSLS']:
        print time, reptiming_box_dict[time][0]
        boxfile.write(reptiming_box_dict[time][0])
        boxfile.write("\n")
    boxfile.close()
    
    os.system("rm -Rf tmp")
    
def ttest(reptiming, rnaseq, gff, max_fpkm, min_fpkm, test_0_expression) :
    print reptiming, rnaseq, gff, max_fpkm, min_fpkm, test_0_expression
    
    min_percentage_overlap = 0
    try :
        os.mkdir("tmp")
    except OSError :
        os.chdir(".")
    
    ### filter FPKM file
    fpkm_list = []
    for a in open(rnaseq, 'r') :
        if a.split("\t")[0] != "tracking_id" :
            fpkm_list.append(float(a.split("\t")[9]))
    filter_file = open("tmp/fpkm.out", 'w')
    for a in open(rnaseq, 'r') :
        if a.split("\t")[0] != "tracking_id" :
            if float(a.split("\t")[9]) > np.percentile(fpkm_list, float(min_fpkm)) and float(a.split("\t")[9]) < np.percentile(fpkm_list, float(max_fpkm)):
                filter_file.write(a)
    filter_file.close()
    
    ### save all FPKM value to a dictionary
    genes = [] 
    gene_fpkm_dict = {}
    random_dict = {}
    for a in open("tmp/fpkm.out", 'r') :
        genes.append(a.split("\t")[0])
        gene_fpkm_dict[a.split("\t")[0]] = float(a.split("\t")[9])
        random_dict[a.split("\t")[0]] = random.uniform(0, 10000)
    
    ### gene coordinates
    genes_bed = open("tmp/gene.bed", 'w')
    gene_coord_dict = {}
    for a in open(gff, 'r') :
        if a[0] != "#" and a.split("\t")[2] == "gene" and a.split("\t")[8].split(";")[0].split(":")[1] in gene_fpkm_dict.keys():
            gene = a.split("\t")[8].split(";")[0].split(":")[1]
            gene_coord_dict[gene] = a.split("\t")[3]+"-"+a.split("\t")[4]
            genes_bed.write(str(a.split("\t")[0]+"\t"+a.split("\t")[3]+"\t"+a.split("\t")[4]+"\t"+gene+"\t"))
            genes_bed.write(str(gene_fpkm_dict[gene])+"\n")
    genes_bed.close()
    
    ### get list of different types of segment
    reptiming_dict = {}
    reptiming_val_dict = {}
    for a in open(reptiming, 'r') :
        if a[0] != "#" :
            reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
            reptiming_val_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
    for a in open(reptiming, 'r') :
        if a[0] != "#" :
            chr = a.split("\t")[0]
            coord_start = int(a.split("\t")[3])-1
            coord_end   = int(a.split("\t")[4])-1
            reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]].append(chr+":"+str(coord_start)+"-"+str(coord_end))
    
    txtfile = open("reptiming_overlap_"+rnaseq.split(".")[0]+".txt", 'w') 
    txtfile.write("Lower FPKM score percentile : "+str(min_fpkm)+"\n")
    txtfile.write("Upper FPKM score percentile : "+str(max_fpkm)+"\n")
    
    ### write into bed and sort
    print "\nstart time:", datetime.datetime.now()
    for time in reptiming_dict.keys() :
        print time
        outfile = open("tmp/"+time+".bed", 'w') 
        for coord in reptiming_dict[time] :
            outfile.write(coord.split(":")[0])
            outfile.write("\t")
            outfile.write(coord.split(":")[1].split("-")[0])
            outfile.write("\t")
            outfile.write(coord.split(":")[1].split("-")[1])
            outfile.write("\t"+str(time)+"\n")
        outfile.close()
        os.system("sortBed -i tmp/"+time+".bed > tmp/"+time+"_sorted.bed")
    
        os.system("intersectBed -wao -a tmp/"+time+".bed -b tmp/gene.bed > tmp/"+time+"_overlap_"+rnaseq.split(".")[0]+".bedGraph")
        overlap_gene = []
        for a in open("tmp/"+time+"_overlap_"+rnaseq.split(".")[0]+".bedGraph",'r') :
            if a.split("\t")[4] != "." :
                overlap = a.split("\t")[-1].strip("\n")
                peak_distance = int(a.split("\t")[2])-int(a.split("\t")[1])
                if float(overlap)/float(peak_distance)*100 >= float(min_percentage_overlap) :
                    overlap_gene.append(a.split("\t")[-3])
    
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
        ax1 = plt.subplot(111)
        plt.hist(t_test_1, bins, alpha=0.5, label='Peak', color='red')
        plt.hist(t_test_0, bins, alpha=0.5, label='No Peak', color='blue')
        plt.xlabel('Log(FPKM)')
        plt.ylabel('Log(frequency)')
        plt.legend(loc='upper right')
        ax1.set_yscale('log')
        ax1.set_yscale('log')
        plt.savefig(time+"_overlap_"+rnaseq.split(".")[0]+"_test.png")
        plt.clf()
        """
        ax2 = plt.subplot(111)
        plt.hist(r_test_1, bins, alpha=0.5, label='Peak', color='red')
        plt.hist(r_test_0, bins, alpha=0.5, label='No Peak', color='blue')
        plt.legend(loc='upper right')
        plt.xlabel('Log(FPKM)')
        plt.ylabel('Log(frequency)')
        ax2.set_yscale('log')
        ax2.set_yscale('log')
        plt.savefig(time+"_overlap_"+rnaseq.split(".")[0]+"_random.png")
        plt.close('all')
        """
    txtfile.write(time+" p-value : "+str(p_val)+" p-value for random : "+str(random_p_val)+"\n")
    txtfile.close()
    os.system("rm -Rf tmp")

parser = OptionParser()
rnagroup = OptionGroup(parser, "Required for ChIP/RNAseq")
rnagroup.add_option("--reptiming", dest="static", help=".gff3 output from Reptiming analysis")
rnagroup.add_option("--rnaseq", dest="rnaseq", help="genes.fpkm_tracking from cufflinks")
rnagroup.add_option("--gff3", dest="gff3", help="Genome annotation gff3 file")
rnagroup.add_option("--max_fpkm", dest="max_fpkm", help="Range maximum percentile of FPKM to filter RNAseq results", default=100)
rnagroup.add_option("--min_fpkm", dest="min_fpkm", help="Range minimum percentile of FPKM to filter RNAseq results", default=0)
rnagroup.add_option("--test_0_expression", action="store_true", dest="test_0_expression", help="test between 0 and > 0 expression",default=False)
rnagroup.add_option("--minoverlap", dest="minoverlap", help="minimum percentage overlap between ChIP/gene (% of peak)", default=0)
rnagroup.add_option("--name", dest="name", help="name for pngs")
parser.add_option_group(rnagroup)
(options, args) = parser.parse_args()

def main():
    #ttest(options.static,options.rnaseq,options.gff3,options.max_fpkm,options.min_fpkm,options.test_0_expression) 
    boxplot(options.static,options.rnaseq,options.gff3,options.max_fpkm,options.min_fpkm,options.minoverlap, options.name)

if __name__ == "__main__":
    main()   
