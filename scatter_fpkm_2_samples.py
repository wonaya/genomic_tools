import os,sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

fpkm_dict = {}

for a in open(sys.argv[1], 'r') :
    if a.split("\t")[0] != "tracking_id" :
        fpkm_dict[a.split("\t")[0]] = []
        fpkm_dict[a.split("\t")[0]].append(float(a.split("\t")[9]))
        
print "finished 1st sample"
        
for b in open(sys.argv[2], 'r') :
    if b.split("\t")[0] != "tracking_id" and b.split("\t")[0] in fpkm_dict.keys() :
        fpkm_dict[b.split("\t")[0]].append(float(b.split("\t")[9]))

star = []
tophat = []
for fpkm in fpkm_dict.keys() :
    if len(fpkm_dict[fpkm]) == 2 :
        tophat.append(fpkm_dict[fpkm][0])
        star.append(fpkm_dict[fpkm][1])
        

plt.scatter(tophat, star, c='red')


plt.xlim(0, 5000)
plt.ylim(0, 5000)



plt.show()

