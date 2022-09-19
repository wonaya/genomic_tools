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

def permutation(static_set, test_set, coordinate_N, permutation_cycle, unique_out, fai, output) :
    try :
        os.mkdir("tmp")
    except OSError :
        os.chdir(".")
    
    os.system("intersectBed -wa -a "+test_set+" -b "+static_set+" > tmp/"+test_set.split(".")[0]+"_overlap_"+static_set.split(".")[0]+".bedGraph")
    a = open("tmp/"+test_set.split(".")[0]+"_overlap_"+static_set.split(".")[0]+".bedGraph",'r')
    alines = a.readlines()
    o = open(test_set, 'r')
    olines = o.readlines()
    original_perc = float(len(set(alines)))/float(len(set(olines)))*100
    print "no. of peaks original:", len(set(alines)),"/", len(set(olines)),":", original_perc
    
    print "\nstart time:", datetime.datetime.now()
    perm_perc = []
    for x in range(0,int(permutation_cycle)) :
        os.system("shuffleBed -i "+test_set+" -g /work/02114/wonaya/genome/maize.genome -excl "+coordinate_N+" > tmp/"+str(x)+".bed")
        os.system("intersectBed -wa -a tmp/"+str(x)+".bed -b "+static_set+" > tmp/"+str(x).split(".")[0]+"_overlap_"+static_set.split(".")[0]+".bedGraph")
        a = open("tmp/"+str(x).split(".bed")[0]+"_overlap_"+static_set.split(".")[0]+".bedGraph", 'r')
        alines = a.readlines()
        perc = float(len(set(alines)))/float(len(set(olines)))*100
        perm_perc.append(perc)
    print "end time:", datetime.datetime.now() # 5 min for 1000 peaks 1000 times, 
    
    count = 0
    for perm in perm_perc :
        if perm >= original_perc : 
            count += 1
    outfile = open(output+".txt", 'w')
    outfile.write("Empirical p-value : "+str(float(count)/float(1000)*100)+"\n")
    outfile.close()
    plt.hist(perm_perc, 50, facecolor='green', alpha=0.75)
    plt.xlabel('percentage overlap (%)')
    plt.ylabel('frequency')
    plt.plot([original_perc, original_perc], [0,100], 'k-', lw=2, color='r')
    plt.savefig(output+".png")
    
    if unique_out != None :
        os.system("intersectBed -v -a "+test_set+" -b "+static_set+" > "+test_set.split(".")[0]+"_unique_from_"+static_set.split(".")[0]+".broadPeak")
        os.system("intersectBed -v -a "+static_set+" -b "+test_set+" > "+static_set.split(".")[0]+"_unique_from_"+test_set.split(".")[0]+".broadPeak")
    os.system("rm -Rf tmp")

parser = OptionParser()
allgroup = OptionGroup(parser, "Required for permutation tests")
allgroup.add_option("--test_1", dest="static", help="not permutated set (.bedGraph)")
allgroup.add_option("--test_2", dest="test", help="permutated set (.bedGraph)")
allgroup.add_option("--Ncoord", dest="ncoord", help="coordinates of N in bedgraph format /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/maize_chrOnly_Ns.bed")
allgroup.add_option("--permutate_n", dest="permutate_n", help="no_of_permutation_cycles_to_run", default=1000)
allgroup.add_option("--unique_out", action="store_false", dest="unique_out", help="output unique peaks as separate bedGraph file")
allgroup.add_option("--fai", dest="fai", help="genome index file")
allgroup.add_option("--output", dest="output", help="output name")
parser.add_option_group(allgroup)
(options, args) = parser.parse_args()

def main():
    permutation(options.static,options.test,options.ncoord,options.permutate_n,options.unique_out,options.fai,options.output) 
              
if __name__ == "__main__":
    main()   
