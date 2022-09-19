import os,sys
import numpy
import multiprocessing
from datetime import datetime
from types import *

def medip_tmp(input, chr, geno):
    ### divide by chromosome
    outfile = open(input+"_chr"+str(chr), 'w')
    for a in open(input, 'r') :
        if a.split("\t")[1] == str(chr) : 
            outfile.write(a)
    outfile.close()
    
    ### divide the file by chr
    print "###start meDIP analysis", chr
    ### meDIP
    chr_dict = {"Mo17":1,"Oh43":2,"CML322":3,"Tx303":4}
    outfile = open(input+"_chr"+str(chr)+"_tmp", 'w')
    for a in open(input+"_chr"+str(chr), 'r') :
        probe_name = []
        b73 = []
        other = []
        for b in open("5geno_fromDiDIP_meDIP_normalized_values_30-1-14.txt" ,'r' ) :
            if b.split("\t")[0] != "chromosome" and b.split("\t")[0] == "chr"+str(chr) :
                if int(b.split("\t")[1]) in range(int(a.split("\t")[2])-300,int(a.split("\t")[3])+301) or int(b.split("\t")[2]) in range(int(a.split("\t")[2])-300,int(a.split("\t")[3])+301) :
                    probe_name.append(b.split("\t")[3])
                    b73.append(float(b.split("\t")[4]))
                    other.append(float(b.split("\t")[4+chr_dict[geno]].strip("\n")))
        if len(probe_name) == 0 :
            outfile.write("\t".join(a.split("\t")[:5]).strip("\r\n")+"\tNA\tNA\tNA\tNA\n")
        else :
            outfile.write("\t".join(a.split("\t")[:5]).strip("\r\n")+"\t"+",".join(probe_name)+"\t"+str(len(probe_name))+"\t"+str(numpy.mean(b73))+"\t"+str(numpy.mean(other))+"\n")
    outfile.close()
"""
### input dmr file, chr, 2nd geno
jobs = []
for chr in range(1,11):
    print chr
    s1 = multiprocessing.Process(target=medip_tmp, args=(sys.argv[1], chr, sys.argv[2], ))
    jobs.append(s1)
    s1.start()
[x.join() for x in jobs]

### merge all

for chr in range(1,11) :
    os.system("cat "+str(sys.argv[1])+"_chr"+str(chr)+"_tmp >> "+str(sys.argv[1])+"_medip_added.txt")
    os.system("rm -Rf "+str(sys.argv[1])+"_chr"+str(chr)+"_tmp")
"""

### calculate each count
outfile = open(str(sys.argv[1]).strip(".txt")+"_medip_added.txt", 'w')
for a in open(str(sys.argv[1])+"_medip_added.txt", 'r'):
    if a.split("\t")[5] != "NA" :
        outfile.write("\t".join(a.split("\t")[:6]))
        outfile.write("\t")
        outfile.write(str(len(a.split("\t")[5].split(","))))
        outfile.write("\t")
        outfile.write("\t".join(a.split("\t")[6:]))
    else :
        outfile.write("\t".join(a.split("\t")[:6]))
        outfile.write("\tNA\tNA\tNA\n")
outfile.close()
    