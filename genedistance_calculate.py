### 

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
#import pandas as pd
from urllib2 import urlopen
import math
import random
from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes

### get list of genes associated with
def permutation(reptiming_gff3, gff, name) :
    print reptiming_gff3, gff, name
    try :
        os.mkdir("tmp")
    except OSError :
        os.chdir(".")
    ### gene coordinates
    genes_bed = open("tmp/gene.bed", 'w')
    gene_coord_dict = {}
    for a in open(gff, 'r') :
        if a[0] != "#" and a.split("\t")[2] == "gene" :
            gene = a.split("\t")[8].split(";")[0].split(":")[1]
            gene_coord_dict[gene] = a.split("\t")[3]+"-"+a.split("\t")[4]
            genes_bed.write(str(a.split("\t")[0]+"\t"+a.split("\t")[3]+"\t"+a.split("\t")[4]+"\t"+gene+"\t"))
            genes_bed.write("0\n")
    genes_bed.close()
    
    ### save gene and gene distance both ways
    os.system("sortBed -i tmp/gene.bed > tmp/gene_sort.bed")
    os.system("mv tmp/gene_sort.bed tmp/gene.bed")
    
    coord_list = []
    gene_list = []
    chr_list = []
    for a in open("tmp/gene.bed", 'r') :
        if a[:2] != "sc" :
            gene_list.append(a.split("\t")[3])
            coord_list.append("-".join(a.split("\t")[:3]))
            chr_list.append(a.split("\t")[0])
    
    gene_distance = {}
    for gene in gene_list :
        gene_distance[gene] = []
    index = 0
    for gene in gene_list :
        #print gene, index, len(gene_list)
        if index == 0 or (coord_list[index-1].split("-")[0]) != (coord_list[index].split("-")[0]) :
            distance_end = int(coord_list[index+1].split("-")[1])-int(coord_list[index].split("-")[2])
            gene_distance[gene].append(distance_end)
            #print distance_end, gene, coord_list[index], coord_list[index+1]
            #sys.exit()
        elif index == len(gene_list)-1 or (coord_list[index+1].split("-")[0]) != (coord_list[index].split("-")[0]) :
            distance_start = int(coord_list[index].split("-")[1])-int(coord_list[index-1].split("-")[2])
            #print distance_start, gene, coord_list[index], coord_list[index-1]
            gene_distance[gene].append(distance_start)
            #sys.exit()
        else :
            distance_start =int(coord_list[index].split("-")[1])-int(coord_list[index-1].split("-")[2])
            distance_end = int(coord_list[index+1].split("-")[1])-int(coord_list[index].split("-")[2])
            gene_distance[gene].extend((distance_start, distance_end))
        #print gene, index, len(gene_list), len(gene_distance[gene])
        
        index += 1
    #print gene_distance
    #sys.exit()
    reptiming_dict = {}
    reptiming_val_dict = {}
    reptiming_count_dict = {}
    for a in open(reptiming_gff3, 'r') :
        if a[0] != "#" :
            reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
            reptiming_val_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
            reptiming_count_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
    for a in open(reptiming_gff3, 'r') :
        if a[0] != "#" :
            chr = a.split("\t")[0]
            coord_start = int(a.split("\t")[3])-1
            coord_end   = int(a.split("\t")[4])-1
            reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]].append(chr+":"+str(coord_start)+"-"+str(coord_end))
    
    ### write Reptiming into bed
    for time in reptiming_dict.keys() :
        #print time
        outfile = open("tmp/"+time+".bed", 'w') 
        for coord in reptiming_dict[time] :
            outfile.write(coord.split(":")[0])
            outfile.write("\t")
            outfile.write(coord.split(":")[1].split("-")[0])
            outfile.write("\t")
            outfile.write(coord.split(":")[1].split("-")[1])
            outfile.write("\t"+str(time)+"\n")
        outfile.close()
        ### intersect reptiming and chip
        os.system("intersectBed -wa -wb -b tmp/gene.bed -a tmp/"+time+".bed > tmp/test_"+time+"_gene_overlap.bedGraph")
        
        sys.exit()
        for a in open("tmp/test_"+time+"_gene_overlap.bedGraph", 'r') :
            if a.split("\t")[-2] in gene_distance.keys() :
                #print a.split("\t")[-2]
                reptiming_count_dict[time].extend(gene_distance[a.split("\t")[-2]])
                print time, a.split("\t")[-2], gene_distance[a.split("\t")[-2]]
        print time, len(reptiming_count_dict[time]), np.mean(reptiming_count_dict[time]), np.percentile(reptiming_count_dict[time], 25), np.median(reptiming_count_dict[time]), np.percentile(reptiming_count_dict[time],75)
        #reptiming_genecount_dict = {}
        #for a in open("tmp/test_"+time+"_gene_overlap.bedGraph", 'r') :
        #    reptiming_genecount_dict["-".join(a.split("\t")[:3])] = 0
        #for a in open("tmp/test_"+time+"_gene_overlap.bedGraph", 'r') :
        #    reptiming_genecount_dict["-".join(a.split("\t")[:3])] += 1
        #total = 0
        #length = len(reptiming_genecount_dict)
        #for dict_keys in reptiming_genecount_dict.keys() :
        #    total += reptiming_genecount_dict[dict_keys]
        #    reptiming_count_dict[time].append(reptiming_genecount_dict[dict_keys])
        #reptiming_val_dict[time] += float(total)/float(length)
    
    
parser = OptionParser()
allgroup = OptionGroup(parser, "Required for permutation tests")
allgroup.add_option("--test_1", dest="reptiming", help="reptiming dataset")
allgroup.add_option("--gff3", dest="gff3", help="gff3 for gene")
allgroup.add_option("--name", dest="name", help="name for plots")
allgroup.add_option("-c", action="store_false", dest="coverage", help="true for calculating coverage instead of count")
parser.add_option_group(allgroup)
(options, args) = parser.parse_args()

def main():
    permutation(options.reptiming,options.gff3,options.name)    
if __name__ == "__main__":
    main()   
 
    
