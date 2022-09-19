import os,sys
import scipy.stats as ss
import numpy as np
import matplotlib.pyplot as plt
import math
import multiprocessing
import readline, glob
import subprocess
from multiprocessing import Process, Manager
from optparse import OptionParser,OptionGroup

"""
parser = OptionParser()
parser.add_option("--runtype", dest="runtype", help="1.bwa, 2.downsample, 3.macs, 4. association")
runtypegrp = OptionGroup(parser, "Required for generating runscripts")
parser.add_option_group(runtypegrp)
runtypegrp.add_option("--input", dest="input", help="Input file for running cutadapt and bwa (space separated) eg.Pair1.fastq Pair2.fastq Barcode")
runtypegrp.add_option("--system", dest="system", help="stampede or lonestar")
runtypegrp.add_option("--time", dest="time", help="no.of hours for bwa (100M reads approx. 4 hours)")
(options, args) = parser.parse_args()

def generater_run(input, system, time):
    ### determine cores required
    a = open(input)
    input_lines = a.readlines()
    cores = len(input_lines)
    runfile = open("temp.txt", 'w') 
    
    outfile = open("run_bwa.sh", 'w')
    if system == "lonestar" : 
        outfile.write("#!/bin/bash\n")
        outfile.write("#$ -V\n")
        outfile.write("#$ -cwd\n")
        outfile.write("#$ -N bwa\n")
        outfile.write("#$ -j y\n")
        outfile.write("#$ -o $JOB_NAME.o$JOB_ID\n")
        outfile.write("#$ -pe 1way "+int(12)*int(cores)+"\n")
        outfile.write("#$ -q normal\n")
        outfile.write("#$ -l h_rt="+str(time)+":00:00\n")
        outfile.write("#$ -A Thompson_Replication\n")
        outfile.write("#$ -P data\n")

def get_dist_plot(name) :
    score_list = []
    for a in open(name+".broadPeak",'r') :
        score_list.append(float(a.split("\t")[7]))

    a = np.array(score_list)
    p5 = np.percentile(a, 5) # return 95th percentile
    p95 = np.percentile(a, 95) # return 95th percentile
    print "5percentile:",p5, "95percentile:",p95
    
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    #ax.set_yscale('log')
    plt.hist(score_list, 100)
    plt.savefig(name+".png")

def s50_peaks(name):    
    outfile = open(name+"_S50.broadPeak",'w')
    for a in open(name+".broadPeak",'r') :
        if float(a.split("\t")[7].strip("\n")) >= 50 :
            outfile.write(a)
    outfile.close()
    
    ### remove redundant lines
    a = open(name+"_S50.broadPeak",'r')
    alines = a.readlines()
    print "no. of peaks", len(set(alines))


def get_gene(chr,return_dict1,return_dict2,return_dict3,return_dict4) :
    distance_set = 1000 ### 1kb

    ### make file to distinguish 5' or 3' direction
    gene_direction_dict = {}
    gene_coord_dict = {}
    for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gff3", 'r') :
        if a.split("\t")[0] == str(chr) and a.split("\t")[2] == "gene" :
            gene_direction_dict[a.split("\t")[8].split(";")[0].split(":")[1]] = a.split("\t")[6]
            gene_coord_dict[a.split("\t")[8].split(";")[0].split(":")[1]] = range(int(a.split("\t")[3]),int(a.split("\t")[4])+1)
            
    ### list of ChIP peaks
    coord_list = []
    for a in open(name+"_S50.broadPeak", 'r') :
        if a.split("\t")[0] == str(chr) :
            coord_list.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
    
    ### how many peaks are within set distance kb of gene
    large_list_of_results = []
    print "total "+name+" peak in chr", chr, ":", len(coord_list)
    
    for coord in coord_list :
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
                elif gene_direction_dict[gene] == "-" : 
                    list_of_results.append("3")
                    list_of_genes.append(gene)
            elif int(coord[0]) > int(gene_coord_dict[gene][-1]) and int(coord[0])-int(gene_coord_dict[gene][-1]) < distance_set :
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
    count_overlap = 0
    count_5 = 0
    count_3 = 0
    count_none = 0
    for data in large_list_of_results :
        count_overlap += data.count('overlap')
        count_5 += data.count('5')
        count_3 += data.count('3')
        count_none += data.count('none')

    print name, "chr"+str(chr), "overlap", count_overlap, "5'", count_5, "3'", count_3, "none", count_none
    
    return_dict1[chr] = count_overlap
    return_dict2[chr] = count_5
    return_dict3[chr] = count_3
    return_dict4[chr] = count_none
    
if options.runtype == "bwa" :
    try :
        generater_run(options.input, options.time, options.system)
    Except IndexError :
        print "hello"
sys.exit()
manager = Manager()
return_dict1 = manager.dict()
return_dict2 = manager.dict()
return_dict3 = manager.dict()
return_dict4 = manager.dict()

jobs = []
for chr in range(1,11) :
    s1 = multiprocessing.Process(target=get_gene, args=(chr,return_dict1, return_dict2, return_dict3, return_dict4))
    jobs.append(s1)
    s1.start()
               
[x.join() for x in jobs]

### ChIP pipeline after running macs2 broadpeak

### length of genome (Zea mays 3.23)
chr_len_dict = {"Pt":140384,"Mt":569630,10:149632204,9:157038028,6:169407836,8:175377492,7:176826311,5:217959525,3:232245527,2:237917468,4:242062272,1:301476924}


b = open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed", 'r') 
blines = b.readlines()
print "no. of genes", len(set(blines))

os.system("intersectBed -wa -a "+name+"_S50.broadPeak -b /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed > "+name+"_Zea_mays.AGPv3.23_protein_coding_exon.bedGraph")

c = open(name+"_Zea_mays.AGPv3.23_protein_coding_exon.bedGraph", 'r')
clines = c.readlines()
print "no. of peaks overlapping with genes", len(set(clines))

"""

### functional association 
### identify unique peaks
def get_2c_8c_unique_peaks(chr,return_dict1,return_dict2) :
    list_2c = []
    for a in open("2C_peaks_S50.broadPeak", 'r') :
        if a.split("\t")[0] == str(chr) :
            list_2c.append(a.split("\t")[1]+"-"+a.split("\t")[2])
    list_8c = []
    for a in open("8C_peaks_S50.broadPeak", 'r') :
        if a.split("\t")[0] == str(chr) :
            list_8c.append(a.split("\t")[1]+"-"+a.split("\t")[2])
    print "starting 2c"
    unique_peaks_2c = []
    for a in list_2c :
        uniq_count = 0
        for b in list_8c : 
            if len(set(range(int(a.split("-")[0]),int(a.split("-")[1])+1))&set(range(int(b.split("-")[0]),int(b.split("-")[1])+1))) > 0 :
                uniq_count += 1
        if uniq_count == 0 :
            unique_peaks_2c.append(a)
    print "starting 8c"
    unique_peaks_8c = []
    for a in list_8c :
        uniq_count = 0
        for b in list_2c : 
            if len(set(range(int(a.split("-")[0]),int(a.split("-")[1])+1))&set(range(int(b.split("-")[0]),int(b.split("-")[1])+1))) > 0 :
                uniq_count += 1
        if uniq_count == 0 :
            unique_peaks_8c.append(a)
    print chr, "2C","8C S50 unique peaks", len(unique_peaks_2c),len(unique_peaks_8c)
    
    return_dict1[chr] = unique_peaks_2c
    return_dict2[chr] = unique_peaks_8c

manager = Manager()
return_dict1 = manager.dict()
return_dict2 = manager.dict()
    
jobs = []
for chr in range(1,11) :
    s1 = multiprocessing.Process(target=get_2c_8c_unique_peaks, args=(chr,return_dict1, return_dict2))
    jobs.append(s1)
    s1.start()

[x.join() for x in jobs]

print len(return_dict1),len(return_dict2)

count = 1
outfile = open("2C_S50_excl_4C_unique_peaks.bedGraph", 'w') 
for chrom in range(1,11) :
    for coord in return_dict1[chrom] :
        outfile.write(str(chrom)+"\t")
        outfile.write(coord.split("-")[0])
        outfile.write("\t")
        outfile.write(coord.split("-")[1])
        outfile.write("\t")
        outfile.write(str(count)+"\n")
        count += 1
outfile.close()
count = 1
outfile = open("8C_S50_excl_4C_unique_peaks.bedGraph", 'w') 
for chrom in range(1,11) :
    for coord in return_dict2[chrom] :
        outfile.write(str(chrom)+"\t")
        outfile.write(coord.split("-")[0])
        outfile.write("\t")
        outfile.write(coord.split("-")[1])
        outfile.write("\t")
        outfile.write(str(count)+"\n")
        count += 1
outfile.close() 
sys.exit()
"""
def get_unique_peaks(chr,return_dict1,return_dict2,return_dict3): 
    list_2c = []
    for a in open("2C_peaks_S50.broadPeak", 'r') :
        if a.split("\t")[0] == str(chr) :
            list_2c.append(a.split("\t")[1]+"-"+a.split("\t")[2])
    list_4c = []
    for a in open("4C_peaks_S50.broadPeak", 'r') :
        if a.split("\t")[0] == str(chr) :
            list_4c.append(a.split("\t")[1]+"-"+a.split("\t")[2])
    list_8c = []
    for a in open("8C_peaks_S50.broadPeak", 'r') :
        if a.split("\t")[0] == str(chr) :
            list_8c.append(a.split("\t")[1]+"-"+a.split("\t")[2])
    print "starting 2c"
    unique_peaks_2c = []
    for a in list_2c :
        uniq_count = 0
        for b in list_4c : 
            if len(set(range(int(a.split("-")[0]),int(a.split("-")[1])+1))&set(range(int(b.split("-")[0]),int(b.split("-")[1])+1))) > 0 :
                uniq_count += 1
        if uniq_count == 0 :
            for c in list_8c : 
                if len(set(range(int(a.split("-")[0]),int(a.split("-")[1])+1))&set(range(int(c.split("-")[0]),int(c.split("-")[1])+1))) > 0 :
                    uniq_count += 1
        if uniq_count == 0 :
            unique_peaks_2c.append(a)
    print "starting 4c"
    unique_peaks_4c = []
    for a in list_4c :
        uniq_count = 0
        for b in list_2c : 
            if len(set(range(int(a.split("-")[0]),int(a.split("-")[1])+1))&set(range(int(b.split("-")[0]),int(b.split("-")[1])+1))) > 0 :
                uniq_count += 1
        if uniq_count == 0 :
            for c in list_8c : 
                if len(set(range(int(a.split("-")[0]),int(a.split("-")[1])+1))&set(range(int(c.split("-")[0]),int(c.split("-")[1])+1))) > 0 :
                    uniq_count += 1
        if uniq_count == 0 :
            unique_peaks_4c.append(a)
    print "starting 8c"
    unique_peaks_8c = []
    for a in list_8c :
        uniq_count = 0
        for b in list_2c : 
            if len(set(range(int(a.split("-")[0]),int(a.split("-")[1])+1))&set(range(int(b.split("-")[0]),int(b.split("-")[1])+1))) > 0 :
                uniq_count += 1
        if uniq_count == 0 :
            for c in list_4c : 
                if len(set(range(int(a.split("-")[0]),int(a.split("-")[1])+1))&set(range(int(c.split("-")[0]),int(c.split("-")[1])+1))) > 0 :
                    uniq_count += 1
        if uniq_count == 0 :
            unique_peaks_8c.append(a)
    print chr, "2C","4C","8C S50 unique peaks", len(unique_peaks_2c),len(unique_peaks_4c),len(unique_peaks_8c)
    
    return_dict1[chr] = unique_peaks_2c
    return_dict2[chr] = unique_peaks_4c
    return_dict3[chr] = unique_peaks_8c

manager = Manager()
return_dict1 = manager.dict()
return_dict2 = manager.dict()
return_dict3 = manager.dict()

    
jobs = []
for chr in range(1,11) :
    s1 = multiprocessing.Process(target=get_unique_peaks, args=(chr,return_dict1, return_dict2, return_dict3))
    jobs.append(s1)
    s1.start()

[x.join() for x in jobs]

print len(return_dict1),len(return_dict2),len(return_dict3)

count = 1
outfile = open("2C_S50_unique_peaks.bedGraph", 'w') 
for chrom in range(1,11) :
    for coord in return_dict1[chrom] :
        outfile.write(str(chrom)+"\t")
        outfile.write(coord.split("-")[0])
        outfile.write("\t")
        outfile.write(coord.split("-")[1])
        outfile.write("\t")
        outfile.write(str(count)+"\n")
        count += 1
outfile.close()
count = 1
outfile = open("4C_S50_unique_peaks.bedGraph", 'w') 
for chrom in range(1,11) :
    for coord in return_dict2[chrom] :
        outfile.write(str(chrom)+"\t")
        outfile.write(coord.split("-")[0])
        outfile.write("\t")
        outfile.write(coord.split("-")[1])
        outfile.write("\t")
        outfile.write(str(count)+"\n")
        count += 1
outfile.close()        
count = 1
outfile = open("8C_S50_unique_peaks.bedGraph", 'w') 
for chrom in range(1,11) :
    for coord in return_dict3[chrom] :
        outfile.write(str(chrom)+"\t")
        outfile.write(coord.split("-")[0])
        outfile.write("\t")
        outfile.write(coord.split("-")[1])
        outfile.write("\t")
        outfile.write(str(count)+"\n")
        count += 1
outfile.close() 

### intersect with gene
os.system("intersectBed -wb -a 2C_S50_unique_peaks.bedGraph -b /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed > 2C_Zea_mays.AGPv3.23_protein_coding_exon.bedGraph")
os.system("intersectBed -wb -a 4C_S50_unique_peaks.bedGraph -b /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed > 4C_Zea_mays.AGPv3.23_protein_coding_exon.bedGraph")
os.system("intersectBed -wb -a 8C_S50_unique_peaks.bedGraph -b /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed > 8C_Zea_mays.AGPv3.23_protein_coding_exon.bedGraph")

gene_list_2c = []
for a in open("2C_Zea_mays.AGPv3.23_protein_coding_exon.bedGraph" ,'r') :
    if a.split("\t")[7].split("_")[0] not in gene_list_2c :
        gene_list_2c.append(a.split("\t")[7].split("_")[0])

outfile = open("gene_list_2c.txt", 'w')
for gene in set(gene_list_2c) :
    outfile.write(gene)
    outfile.write("\n")
outfile.close()

gene_list_4c = []
for a in open("4C_Zea_mays.AGPv3.23_protein_coding_exon.bedGraph" ,'r') :
    if a.split("\t")[7].split("_")[0] not in gene_list_2c :
        gene_list_4c.append(a.split("\t")[7].split("_")[0])

outfile = open("gene_list_4c.txt", 'w')
for gene in set(gene_list_4c) :
    outfile.write(gene)
    outfile.write("\n")
outfile.close()

gene_list_8c = []
for a in open("8C_Zea_mays.AGPv3.23_protein_coding_exon.bedGraph" ,'r') :
    if a.split("\t")[7].split("_")[0] not in gene_list_2c :
        gene_list_8c.append(a.split("\t")[7].split("_")[0])

outfile = open("gene_list_8c.txt", 'w')
for gene in set(gene_list_8c) :
    outfile.write(gene)
    outfile.write("\n")
outfile.close()


 
pop_count = len(set(pop_list))
### make assoc file Gene:Go1;Go2
if not os.path.isfile("gene_goterm.txt") :
    gene_dict = {}
    a = open("/work/02114/wonaya/genome/annotation/MapMan_B73_annot.txt", 'r')
    lines = a.readlines()
    for line in lines[1:] :
        gene_dict[line.split("\t")[0]] = []
    print len(gene_dict)
    for line in lines[1:] :
        if line.split("\t")[1].strip("\n") not in gene_dict[line.split("\t")[0]] :
            gene_dict[line.split("\t")[0]].append(line.split("\t")[1].strip("\n"))
    outfile = open("gene_goterm.txt", 'w') 
    for gene in gene_dict.keys() :
        outfile.write(gene+"\t"+";".join(gene_dict[gene])+"\n")
    outfile.close()

### make description file
if not os.path.isfile("go_desc.txt") :
    go_dict = {}
    a = open("/work/02114/wonaya/genome/annotation/Zm_B73_5b_FGS_cds_2012_mapman_summary.txt", 'r')
    lines = a.readlines()
    for line in lines[1:] :
        go_dict[line.split("\t")[0]] = line.split("\t")[1]
    print len(go_dict)
    outfile = open("go_desc.txt", 'w') 
    for go in go_dict.keys() :
        outfile.write(go+"\t"+go_dict[go]+"\n")
    outfile.close()
"""

