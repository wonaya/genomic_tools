### TO DO: UPDATE NAME SO THAT OVERLAP DOES NOT CONFUSE OTHERS
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

def permutation(static_set, test_set, coordinate_N, fai, permutation_cycle, max_macs, min_macs) :
    print static_set, test_set, coordinate_N, fai, permutation_cycle, max_macs, min_macs
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
    ### second ChIP
    chip_val_list = []
    for a in open(static_set, 'r') :
        chip_val_list.append(float(a.split("\t")[6]))
    filter_file = open("tmp/static_chip.bed", 'w')
    for a in open(static_set, 'r') :
        if float(a.split("\t")[6]) > np.percentile(chip_val_list, float(min_macs)) and float(a.split("\t")[6]) < np.percentile(chip_val_list, float(max_macs)):
            filter_file.write(a)
    filter_file.close()
    
    #outfile = open(test_set.split(".")[0]+"_overlap_"+static_set.split(".")[0]+".txt", 'w')
    #outfile.write("Lower MACS score percentile : "+str(min_macs)+"\n")
    #outfile.write("Upper MACS score percentile : "+str(max_macs)+"\n")
    #outfile.write("No. of permutation cycles : "+str(permutation_cycle)+"\n")
    
    ### intersect 
    os.system("intersectBed -wa -a tmp/test_chip.bed -b tmp/static_chip.bed > tmp/test_static_overlap.bedGraph")
    a = open("tmp/test_static_overlap.bedGraph",'r')
    alines = a.readlines()
    o = open("tmp/test_chip.bed", 'r')
    olines = o.readlines()
    original_perc = float(len(set(alines)))/float(len(set(olines)))*100
    outfile.write("No. of peaks original: "+str(len(set(alines)))+"/"+str(len(set(olines)))+" : "+str(original_perc)+"\n")
    
    print "\nstart time:", datetime.datetime.now()
    perm_perc = []
    for x in range(0,int(permutation_cycle)) :
        os.system("shuffleBed -i tmp/test_chip.bed -g "+fai+" -excl "+coordinate_N+" > tmp/"+str(x)+".bed")
        os.system("intersectBed -wa -a tmp/"+str(x)+".bed -b tmp/static_chip.bed > tmp/"+str(x).split(".")[0]+"_overlap_static.bedGraph")
        a = open("tmp/"+str(x).split(".")[0]+"_overlap_static.bedGraph", 'r')
        alines = a.readlines()
        perc = float(len(set(alines)))/float(len(set(olines)))*100
        perm_perc.append(perc)
    print "end time:", datetime.datetime.now() # 5 min for 1000 peaks 1000 times, 
    
    count = 0
    for perm in perm_perc :
        if perm >= original_perc : 
            count += 1
    
    #outfile.write("Empirical p-value : "+str(float(count)/float(1000)*100)+"\n")
    #outfile.close()
    plt.hist(perm_perc, 50, facecolor='green', alpha=0.75)
    plt.xlabel('percentage overlap (%)')
    plt.ylabel('frequency')
    plt.plot([original_perc, original_perc], [0,100], 'k-', lw=2, color='r')
    
    plt.savefig(test_set.split(".")[0]+"_overlap_"+static_set.split(".")[0]+".png")
    
    ### unique from bedgraph
    os.system("intersectBed -v -a tmp/test_chip.bed -b tmp/static_chip.bed > "+test_set.split(".")[0]+"_unique_from_"+static_set.split(".")[0]+".broadPeak")
    os.system("intersectBed -v -a tmp/static_chip.bed -b tmp/test_chip.bed > "+static_set.split(".")[0]+"_unique_from_"+test_set.split(".")[0]+".broadPeak")
    
    ### overlap bedGraph
    os.system("intersectBed -wa -a tmp/test_chip.bed -b tmp/static_chip.bed > "+test_set.split(".")[0]+"_overlap_"+static_set.split(".")[0]+".broadPeak")
    
    os.system("rm -Rf tmp")

parser = OptionParser()
allgroup = OptionGroup(parser, "Required for permutation tests")
allgroup.add_option("--test_1", dest="static", help="not permutated set (.broadPeak)")
allgroup.add_option("--test_2", dest="test", help="permutated set (.broadPeak)")
allgroup.add_option("--Ncoord", dest="ncoord", help="coordinates of N in bedgraph format maize_chrOnly_Ns.bed")
allgroup.add_option("--permutate_n", dest="permutate_n", help="no_of_permutation_cycles_to_run", default=1000)
allgroup.add_option("--fai", dest="fai", help="genome index file")
allgroup.add_option("--max_macs", dest="max_macs", help="Range maximum percentile of macs scores to filter MACS scores for both ChIP")
allgroup.add_option("--min_macs", dest="min_macs", help="Range minimum percentile of macs scores to filter MACS scores for both ChIP")
parser.add_option_group(allgroup)
(options, args) = parser.parse_args()

def main():
    permutation(options.static,options.test,options.ncoord,options.fai,options.permutate_n,options.max_macs,options.min_macs)    
if __name__ == "__main__":
    main()   
