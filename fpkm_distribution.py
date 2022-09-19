import os,sys
import scipy.stats as ss
import numpy as np
import matplotlib.pyplot as plt
import math
### usage fpkm_distribution.py genes.fpkm_tracking


## fpkm data range
fpkm_list = []
for a in open(sys.argv[1],'r') :
    if a.split("\t")[0] != "tracking_id" and float(a.split("\t")[9]) > 0 :
        gene = a.split("\t")[0]
        #fpkm = math.log(float(a.split("\t")[9]),10)
        fpkm_list.append(float(a.split("\t")[9]))

a = np.array(fpkm_list)
p25 = np.percentile(a, 25) # return 95th percentile
p50 = np.percentile(a, 50) # return 95th percentile
p75 = np.percentile(a, 75) # return 95th percentile
p100 = np.percentile(a, 100) # return 95th percentile
print p25, p50, p75, p100
#sys.exit()

fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.set_yscale('log')
plt.hist(fpkm_list, 100)
plt.show()    
"""
## fpkm to bedGraph

outfile = open("genes.fpkm.bedGraph", 'w') 
for a in open(sys.argv[1], 'r') :
    if a.split("\t")[0] != "tracking_id" and float(a.split("\t")[9]) > 0 :
        outfile.write(a.split("\t")[6].split(":")[0])
        outfile.write("\t")
        outfile.write(a.split("\t")[6].split(":")[1].split("-")[0])
        outfile.write("\t")
        outfile.write(a.split("\t")[6].split(":")[1].split("-")[1])
        outfile.write("\t")
        outfile.write(a.split("\t")[9])
        outfile.write("\t")
        outfile.write(a.split("\t")[0])
        outfile.write("\n")
outfile.close()
        
"""
