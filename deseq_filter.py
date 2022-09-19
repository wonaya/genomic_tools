import os,sys

outfile = open("valid_embryo_vs_root.txt", 'w')
pval_list = []
for a in open("embryo_vs_root.txt", 'r') :
    if a.split("\t")[2] != '"baseMeanA"' :
        if float(a.split("\t")[3]) != 0 and float(a.split("\t")[4]) != 0:
            pval_list.append(float(a.split("\t")[7]))
            outfile.write(a)
outfile.close()

### draw histogram of pvalues
import numpy as np
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
numBins = 50
ax.hist(pval_list,numBins,color='green',alpha=0.8)
#plt.show()

outfile = open("valid_sig_embryo_vs_root.txt", 'w')
pval_sig_list = []
log_fold_list = []
for a in open("valid_embryo_vs_root.txt", 'r') :
    if float(a.split("\t")[7]) <= 1e-20 :
        pval_sig_list.append(float(a.split("\t")[7]))
        log_fold_list.append(float(a.split("\t")[6]))
        outfile.write(a)
outfile.close()

fig = plt.figure()
ax = fig.add_subplot(111)
numBins = 50
ax.hist(pval_sig_list,numBins,color='green',alpha=0.8)
#plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
numBins = 50
ax.hist(log_fold_list,numBins,color='green',alpha=0.8)
#plt.show()

outfile = open("valid_sig_upinroot_vs_embryo.txt", 'w')
for a in open("valid_sig_embryo_vs_root.txt", 'r') :
    if float(a.split("\t")[6]) > 5 :
        outfile.write(a)
outfile.close()

outfile = open("valid_sig_upinembryo_vs_root.txt", 'w')
for a in open("valid_sig_embryo_vs_root.txt", 'r') :
    if float(a.split("\t")[6]) < -5 :
        outfile.write(a)
outfile.close()

### extract gene list only
outfile = open("valid_sig_upinembryo_vs_root_gene_list.txt",'w')
for a in open("valid_sig_upinembryo_vs_root.txt",'r') :
    outfile.write(a.split("\t")[1].strip('"'))
    outfile.write("\n")
outfile.close() 
outfile = open("valid_sig_upinroot_vs_embryo_gene_list.txt",'w')
for a in open("valid_sig_upinroot_vs_embryo.txt",'r') :
    outfile.write(a.split("\t")[1].strip('"'))
    outfile.write("\n")
outfile.close()       

### run gsea.py