import os,sys
import math

## nearest CG, CHG
list_loc = []
for a in open("CpG_context_B73_all3_R1_bt202_chr10.bismark.cov", 'r') :
    list_loc.append(int(a.split("\t")[1]))
        
list_distance = []
index = 0
for loc in list_loc :
    if index > 0 and index < len(list_loc)-1 :
        small_list = []
        small_list.append(loc-list_loc[index-1])
        small_list.append(list_loc[index+1]-loc)
        list_distance.append(math.log(min(small_list),2))
    elif index == 0 :
        list_distance.append(math.log(list_loc[index+1]-loc,2))
    elif index == len(list_loc)-1 :
        list_distance.append(math.log(loc-list_loc[index-2],2))
    index += 1

import matplotlib.pyplot as plt

#xbins = [0, len(list_distance)]
#plt.hist(x, bins=xbins, color = 'blue') 
#Does not make the histogram correct. It counts the occurances of the individual counts. 
fig, ax = plt.subplots()
ax.set_yscale('log',basey=2)
n, bins, patches = ax.hist(list_distance, 50, facecolor='green', alpha=0.75)
plt.xlabel('Log2 Distance')
plt.ylabel('Log2 Frequency')

plt.yscale('log')
#plot works but I need this in histogram format
#plt.show()
plt.savefig("CpG_B73_chr10_log2distance_log2freq.png")