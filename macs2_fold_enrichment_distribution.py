### distribution of macs2 scores (fold enrichment)

import os,sys
import matplotlib.pyplot as plt
import numpy as np

macs_fold_list = []
logp_fold_list = []
logq_fold_list = []
for a in open(sys.argv[1], 'r') : ### .broadPeak
    macs_fold_list.append(float(a.split("\t")[-3]))
    logp_fold_list.append(float(a.split("\t")[-2]))
    logq_fold_list.append(float(a.split("\t")[-1].strip("\n")))
bins = np.linspace(0,100, 100) #min, max
    
plt.hist(macs_fold_list, bins, alpha=0.5, label='Fold')
plt.hist(logp_fold_list, bins, alpha=0.5, label='-log(P)')
plt.hist(logq_fold_list, bins, alpha=0.5, label='-log(Q)')
plt.legend(loc='upper right')
    
plt.show()
