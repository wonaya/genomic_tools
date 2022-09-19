import os,sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import matplotlib
import math
matplotlib.use('Agg')

#### FPKM distribution

fpkm_list = []
count_0 = 0
for a in open(sys.argv[1], 'r') :
    if a.split("\t")[0] != "tracking_id" :
        if float(a.split("\t")[9]) == 0 :
            #count_0 += 1
            fpkm_list.append(float(a.split("\t")[9]))
        else :
            #fpkm_list.append(math.log(float(a.split("\t")[9]),2))
            fpkm_list.append(float(a.split("\t")[9]))
print "count 0", count_0
    
bins = np.linspace(0, max(fpkm_list), 100)

ax1 = plt.subplot(111)
plt.hist(fpkm_list, bins, alpha=0.5, label='Peak', color='red')
plt.legend(loc='upper right')
plt.xlabel('Log(FPKM)')
plt.ylabel('Log(frequency)')
ax1.set_yscale('log')
plt.show()
plt.close('all')
