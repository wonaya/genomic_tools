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

def peak_distribution(name, score) :
    score_list = []
    for a in open(name,'r') :
        score_list.append(float(a.split("\t")[7]))

    a = np.array(score_list)
    p5 = np.percentile(a, 5) # return 95th percentile
    p95 = np.percentile(a, 95) # return 95th percentile
    print name,":   5 percentile:",p5, "95 percentile:",p95
    
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    #ax.set_yscale('log')
    plt.hist(score_list, 100)
    plt.savefig(name.split(".")[0]+".png")
    
    outfile = open(name.split(".")[0]+"_S"+str(score)+".broadPeak",'w')
    for a in open(name,'r') :
        if float(a.split("\t")[7].strip("\n")) >= score :
            outfile.write(a)
    outfile.close()
    
    a = open(name.split(".")[0]+"_S"+str(score)+".broadPeak",'r')
    alines = a.readlines()
    print name, score, "no. of peaks", len(set(alines))

def unique_peaks(list_of_bams, fai) :
    cores = float(multiprocessing.cpu_count())
    chr_list = []
    for a in open(fai, 'r') :
        chr_list.append(a.split("\t")[0])
    total_rounds = len(chr_list)/cores
    
    manager = Manager()
    dict_list = []
    return_dict = manager.dict()
    
    for rounds in range(0, int(total_rounds)) :
        jobs = []
        result_queue = multiprocessing.Queue()
        for chr in chr_list[int(int(cores)*rounds):int(int(cores)*rounds+cores)] :
        #for chr in ["10"]:
            s1 = multiprocessing.Process(target=get_unique_peaks, args=(list_of_bams,chr,return_dict))
            jobs.append(s1)
            s1.start()
        [x.join() for x in jobs]
    
def get_unique_peaks(list_of_bams,chr,return_dict) :
    large_list = []
    for bam in list_of_bams.split(",") :
        list_test = []
        for a in open(bam, 'r') :
            if a.split("\t")[0] == str(chr) :
                list_test.append(a.split("\t")[1]+"-"+a.split("\t")[2])
        large_list.append(list_test)
    unique_peaks_test_list = []
    for bam in list_of_bams.split(",") :
        unique_peaks_test = []
        for list in large_list :
            if list_of_bams.split(",").index(bam) != large_list.index(list) :
                for a in large_list[list_of_bams.split(",").index(bam)] :
                    uniq_count = 0
                    for b in list : 
                        if len(set(range(int(a.split("-")[0]),int(a.split("-")[1])+1))&set(range(int(b.split("-")[0]),int(b.split("-")[1])+1))) > 0 :
                            uniq_count += 1
                    if uniq_count == 0 :
                        unique_peaks_test.append(a)
        unique_peaks_test_list.append(unique_peaks_test)
    return_dict[chr] = unique_peaks_test_list
    return return_dict
                
parser = OptionParser()
parser.add_option("--runtype", dest="runtype", help="1.filter 2.unique")
peakfilter = OptionGroup(parser, "Required for running filtering peaks")
parser.add_option_group(peakfilter)
peakfilter.add_option("--broadpeak", dest="broadpeak", help="broadpeak file to test")
peakfilter.add_option("--score", dest="score", help="macs2 score to filter from (default=50)", default=50)
findunique = OptionGroup(parser, "Required for identifying unique peaks")
parser.add_option_group(findunique)
findunique.add_option("--listofpeaks", dest="listofpeaks", help="list of peak files separated by comma(,)")
findunique.add_option("--fai", dest="fai", help="fasta index file")
(options, args) = parser.parse_args()

if options.runtype == "filter" :
    peak_distribution(options.broadpeak,options.score)

elif options.runtype == "unique" :
    unique_peaks(options.listofpeaks,options.fai)

