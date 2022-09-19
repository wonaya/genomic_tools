### repk4k27compare.py

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

def chipchiptimecompare(reptiming_gff3, test_set, test2_set, name) :
    print reptiming_gff3, test_set, test2_set, name
    ### create tmp directory
    try :
        os.mkdir("tmp")
    except OSError :
        os.chdir(".")
    
    ### initial information written into file
    txtfile = open(name+"_overlap_reptiming.txt", 'w')
    
    ### filter set to min and max range
    ### first ChIP K4
    filter_file = open("tmp/test1_chip.bed", 'w')
    for a in open(test_set, 'r') :
        filter_file.write(a)
    filter_file.close()
    
    ### second ChIP K27
    filter_file = open("tmp/test2_chip.bed", 'w')
    for a in open(test2_set, 'r') :
        filter_file.write(a)
    filter_file.close()
        
    ### prepare dictionaries for reptiming
    reptiming_dict = {}
    reptiming_val_dict = {}
    reptiming_count_dict = {}
    reptiming_total_dict = {}
    reptiming_chip1_dict = {}
    reptiming_chip2_dict = {}
    for a in open(reptiming_gff3, 'r') :
        if a[0] != "#" :
            reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
            reptiming_val_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
            reptiming_count_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
            reptiming_total_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
            reptiming_chip1_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
            reptiming_chip2_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = 0
    
    
    for a in open(reptiming_gff3, 'r') :
        if a[0] != "#" :
            chr = a.split("\t")[0]
            coord_start = int(a.split("\t")[3])-1
            coord_end   = int(a.split("\t")[4])-1
            reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]].append(chr+":"+str(coord_start)+"-"+str(coord_end))
    
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
        
        ### overlap first chip
        os.system("intersectBed -wa -b tmp/test1_chip.bed -a tmp/"+time+".bed > tmp/test1_"+time+"_overlap.bedGraph")
        
        ### overlap second chip
        os.system("intersectBed -wa -b tmp/test2_chip.bed -a tmp/"+time+".bed > tmp/test2_"+time+"_overlap.bedGraph")
        
        ### get regions where there's both chip1 and chip2 peaks in each segment
        afile = open("tmp/test1_"+time+"_overlap.bedGraph", 'r') 
        alines = afile.readlines()
        bfile = open("tmp/test2_"+time+"_overlap.bedGraph", 'r') 
        blines = bfile.readlines()
        
        count = 0
        total = 0
        for a in open("tmp/"+time+".bed", 'r') :
            total += 1
            if a in alines and a in blines :
                count += 1
        chip1_count = 0
        for a in open("tmp/"+time+".bed", 'r') :
            if a in alines :
                chip1_count += 1
        chip2_count = 0
        for a in open("tmp/"+time+".bed", 'r') :
            if a in blines :
                chip2_count += 1
            
        ### both
        reptiming_count_dict[time] += count
        ### total
        reptiming_total_dict[time] += total
        ### chip1
        reptiming_chip1_dict[time] += chip1_count
        ### chip2
        reptiming_chip2_dict[time] += chip2_count
            
    count_val = [reptiming_count_dict['ES'],reptiming_count_dict['ESMS'],reptiming_count_dict['MS'],reptiming_count_dict['MSLS'],reptiming_count_dict['LS'],reptiming_count_dict['ESLS'],reptiming_count_dict['ESMSLS']]
    total_val = [reptiming_total_dict['ES'],reptiming_total_dict['ESMS'],reptiming_total_dict['MS'],reptiming_total_dict['MSLS'],reptiming_total_dict['LS'],reptiming_total_dict['ESLS'],reptiming_total_dict['ESMSLS']]
    chip1_val = [reptiming_chip1_dict['ES'],reptiming_chip1_dict['ESMS'],reptiming_chip1_dict['MS'],reptiming_chip1_dict['MSLS'],reptiming_chip1_dict['LS'],reptiming_chip1_dict['ESLS'],reptiming_chip1_dict['ESMSLS']]
    chip2_val = [reptiming_chip2_dict['ES'],reptiming_chip2_dict['ESMS'],reptiming_chip2_dict['MS'],reptiming_chip2_dict['MSLS'],reptiming_chip2_dict['LS'],reptiming_chip2_dict['ESLS'],reptiming_chip2_dict['ESMSLS']]
    
    print count_val
    print total_val
    print chip1_val
    print chip2_val
    
    ind = np.arange(len(reptiming_count_dict))   
    width= 0.25
    fig,  ax= plt.subplots()
    rects1 = ax.bar(ind, count_val, width,color='r')
    rects2 = ax.bar(ind+width, chip1_val, width,color='b')
    rects3 = ax.bar(ind+width+width, chip2_val, width,color='g')
    rects4 = ax.bar(ind+width+width+width, total_val, width,color='y')
    ax.set_ylabel('Segment Count')
    ax.set_xticks(ind+width+width)
    ax.set_xticklabels(('ES','ESMS','MS','MSLS','LS','ESLS','ESMSLS'))
    ax.legend((rects1[0], rects2[0],rects3[0], rects4[0]), ('Segments with K4me3,K27me3', 'Segments with K4me3','Segments with K27me3','Total segments'))
    
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,'%d' % int(height),ha='center', va='bottom')
    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)
    autolabel(rects4)
    #plt.figure(figsize=(20,10))
    plt.savefig(name+'.png')
    plt.close()
    
parser = OptionParser()
allgroup = OptionGroup(parser, "Required for permutation tests")
allgroup.add_option("--test_1", dest="reptiming", help="reptiming dataset")
allgroup.add_option("--test_2", dest="chip1", help="1st chip peaks (.broadPeak)")
allgroup.add_option("--test_3", dest="chip2", help="2nd chip peaks (.broadPeak)")
allgroup.add_option("--name", dest="name", help="Name for the plot")
parser.add_option_group(allgroup)
(options, args) = parser.parse_args()

def main():
    chipchiptimecompare(options.reptiming,options.chip1,options.chip2,options.name)    

if __name__ == "__main__":
    main()   
 
