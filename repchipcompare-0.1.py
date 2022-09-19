### repchipcompare-0.1.py

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

def permutation(reptiming_gff3, test_set, coordinate_N, fai, permutation_cycle, max_macs, min_macs, name) :
    print reptiming_gff3, test_set, coordinate_N, fai, permutation_cycle, max_macs, min_macs, name
    try :
        os.mkdir("tmp")
    except OSError :
        os.chdir(".")
    ### filter set to min and max range
    ### first ChIP
    chip_val_list = []
    for a in open(test_set, 'r') :
        chip_val_list.append(float(a.split("\t")[6]))
    filter_file = open("tmp/test_chip.bed", 'w')
    for a in open(test_set, 'r') :
        if float(a.split("\t")[6]) > np.percentile(chip_val_list, float(min_macs)) and float(a.split("\t")[6]) < np.percentile(chip_val_list, float(max_macs)):
            filter_file.write(a)
    filter_file.close()
    
    ### initial information written into file
    txtfile = open(test_set.split(".")[0]+"_overlap_reptiming.txt", 'w')
    txtfile.write("Lower MACS score percentile : "+str(min_macs)+"\n")
    txtfile.write("Upper MACS score percentile : "+str(max_macs)+"\n")
    txtfile.write("No. of permutation cycles : "+str(permutation_cycle)+"\n")
    
    ### get list of different types of segment
    reptiming_dict = {}
    reptiming_val_dict = {}
    reptiming_count_dict = {}
    reptiming_macs_dict = {}
    reptiming_box_dict = {}
    for a in open(reptiming_gff3, 'r') :
        if a[0] != "#" :
            reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
            reptiming_val_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
            reptiming_count_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
            reptiming_macs_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
            reptiming_box_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
    for a in open(reptiming_gff3, 'r') :
        if a[0] != "#" :
            chr = a.split("\t")[0]
            coord_start = int(a.split("\t")[3])-1
            coord_end   = int(a.split("\t")[4])-1
            reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]].append(chr+":"+str(coord_start)+"-"+str(coord_end))
    ### write into bed and sort
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
        os.system("sortBed -i tmp/"+time+".bed > tmp/"+time+"_sorted.bed")
    
        ### intersect 
        os.system("intersectBed -wa -a tmp/test_chip.bed -b tmp/"+time+".bed > tmp/test_"+time+"_overlap.bedGraph")
        for a in open("tmp/test_"+time+"_overlap.bedGraph", 'r') :
            reptiming_macs_dict[time].append(float(a.split("\t")[-2]))
        B= plt.boxplot(reptiming_macs_dict[time])
        reptiming_box_dict[time].append(str(time)+"\t"+str([item.get_ydata()[1] for item in B['whiskers']][0])+"\t"+str([item.get_ydata()[0] for item in B['whiskers']][0])+"\t"+str(np.median(reptiming_macs_dict[time]))+"\t"+str([item.get_ydata()[0] for item in B['whiskers']][1])+"\t"+str([item.get_ydata()[1] for item in B['whiskers']][1]))
        a = open("tmp/test_"+time+"_overlap.bedGraph",'r')
        alines = a.readlines()
        o = open("tmp/test_chip.bed", 'r')
        olines = o.readlines()
        original_perc = float(len(set(alines)))/float(len(set(olines)))*100
        reptiming_count_dict[time] += len(set(alines))
    
        txtfile.write(time+": no. of peaks original: "+str(len(set(alines)))+"/"+str(len(set(olines)))+" : "+str(original_perc)+"\n")
        reptiming_val_dict[time] = original_perc
        perm_perc = []
        for x in range(0,int(permutation_cycle)) :
            os.system("shuffleBed -i tmp/test_chip.bed -g "+fai+" -excl "+coordinate_N+" > tmp/"+str(x)+".bed")
            os.system("intersectBed -wa -a tmp/"+str(x)+".bed -b tmp/"+time+".bed > tmp/"+str(x).split(".")[0]+"_overlap_static.bedGraph")
            a = open("tmp/"+str(x).split(".")[0]+"_overlap_static.bedGraph", 'r')
            alines = a.readlines()
            perc = float(len(set(alines)))/float(len(set(olines)))*100
            perm_perc.append(perc)
        count = 0
        for perm in perm_perc :
            if perm >= original_perc : 
                count += 1
        txtfile.write(time+": empirical p-value : "+str((float(count)/float(permutation_cycle)))+"\n")
        
        plt.hist(perm_perc, 50, facecolor='green', alpha=0.75)
        plt.xlabel('percentage overlap (%)')
        plt.ylabel('frequency')
        plt.plot([original_perc, original_perc], [0,100], 'k-', lw=2, color='r')
        plt.savefig(test_set.split(".")[0]+"_overlap_"+time+".png")
        plt.close('all')
    plt.bar(range(len(reptiming_val_dict)), reptiming_val_dict.values(), align='center')
    plt.xticks(range(len(reptiming_val_dict)), reptiming_val_dict.keys())
    plt.show()
    txtfile.close()
    
    ### bar graph showing ChIP peaks count per segment
    val = [reptiming_count_dict['ES'],reptiming_count_dict['ESMS'],reptiming_count_dict['MS'],reptiming_count_dict['MSLS'],reptiming_count_dict['LS'],reptiming_count_dict['ESLS'],reptiming_count_dict['ESMSLS']]
    
    ind = np.arange(len(reptiming_count_dict))   
    width= 0.35
    fig,  ax= plt.subplots()
    rects1 = ax.bar(ind, val, width,color='r')
    ax.set_ylabel('Peaks count')
    #ax.set_xticks(ind+width)
    ax.set_xticklabels(('ES','ESMS','MS','MSLS','LS','ESLS','ESMSLS'))
    #ax.legend((rects1[0], rects2[0]), ('>0', '=0'))
    
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,'%d' % int(height),ha='center', va='bottom')
    autolabel(rects1)
    plt.savefig(name+'_barchart.png')
    plt.close()
    
    es_test = np.array(reptiming_macs_dict['ES'])
    esms_test = np.array(reptiming_macs_dict['ESMS'])
    ms_test = np.array(reptiming_macs_dict['MS'])
    msls_test = np.array(reptiming_macs_dict['MSLS'])
    ls_test = np.array(reptiming_macs_dict['LS'])
    esls_test = np.array(reptiming_macs_dict['ESLS'])
    esmsls_test = np.array(reptiming_macs_dict['ESMSLS'])
    data = [es_test,esms_test,ms_test,msls_test,ls_test,esls_test,esmsls_test]
    fig = plt.figure()
    ax1= fig.add_subplot(1,1,1)
    ax1.boxplot(data,showfliers=False)
    ax1.set_ylabel('MACS score')
    ax1.set_xticklabels(['ES', 'ESMS','MS', 'MSLS', 'LS','ESLS', 'ESMSLS'])
    ax1.get_xaxis().tick_bottom()
    plt.savefig(name+'_boxplot.png')
    outfile = open(name+"_boxplot_value.txt", 'w')
    for time in ['ES', 'ESMS','MS', 'MSLS', 'LS','ESLS', 'ESMSLS']:
        print reptiming_box_dict[time][0]
        outfile.write(reptiming_box_dict[time][0])
        outfile.write("\n")
    outfile.close()
     ### unique from bedgraph
    #os.system("intersectBed -v -a tmp/test_chip.bed -b tmp/static_chip.bed > "+test_set.split(".")[0]+"_unique_from_"+static_set.split(".")[0]+".broadPeak")
    #os.system("intersectBed -v -a tmp/static_chip.bed -b tmp/test_chip.bed > "+static_set.split(".")[0]+"_unique_from_"+test_set.split(".")[0]+".broadPeak")
    
    ### overlap bedGraph
    #os.system("intersectBed -wa -a tmp/test_chip.bed -b tmp/static_chip.bed > "+test_set.split(".")[0]+"_overlap_"+static_set.split(".")[0]+".broadPeak")
    
    os.system("rm -Rf tmp")

parser = OptionParser()
allgroup = OptionGroup(parser, "Required for permutation tests")
allgroup.add_option("--test_1", dest="reptiming", help="reptiming dataset")
allgroup.add_option("--test_2", dest="chip", help="chip peaks (.broadPeak)")
allgroup.add_option("--Ncoord", dest="ncoord", help="coordinates of N in bedgraph format maize_chrOnly_Ns.bed")
allgroup.add_option("--permutate_n", dest="permutate_n", help="no_of_permutation_cycles_to_run", default=10)
allgroup.add_option("--fai", dest="fai", help="genome index file")
allgroup.add_option("--max_macs", dest="max_macs", default=100, help="Range maximum percentile of macs scores to filter MACS scores for both ChIP")
allgroup.add_option("--min_macs", dest="min_macs", default=0, help="Range minimum percentile of macs scores to filter MACS scores for both ChIP")
allgroup.add_option("--name", dest="name", help="name for plots")
parser.add_option_group(allgroup)
(options, args) = parser.parse_args()

def main():
    permutation(options.reptiming,options.chip,options.ncoord,options.fai,options.permutate_n,options.max_macs,options.min_macs, options.name)    
if __name__ == "__main__":
    main()   
