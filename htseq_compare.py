import os,sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import * 

## use htseq_compare.py control.htseq.out file1 file2 ... fileN

nooffiles = len(sys.argv)

control = []
a0 = open(sys.argv[1], 'r')
a0lines = a0.readlines()

for a0line in a0lines :
    control.append(int(a0line.split("\t")[1].strip("\n")))

large_list = []
for x in range(2,nooffiles) :
#for x in range(2,3):
    sample = []
    a1 = open(sys.argv[x], 'r')
    a1lines = a1.readlines()
    for a1line in a1lines :
        sample.append(int(a1line.split("\t")[1].strip("\n")))
    large_list.append(sample)

print len(control), len(large_list[0]),len(large_list[1]),len(large_list[2]),len(large_list[3])
fig = plt.figure(figsize=(8, 8))
a1 = fig.add_subplot(2, 2, 1)
a1.set_xscale('log')
a1.set_yscale('log')
m,b = polyfit(control[:-5], large_list[0][:-5], 1) 
plt.scatter(control[:-5], large_list[0][:-5],color='blue')
plt.plot(range(1,10000),range(1,10000), 'k-')
plt.xlim(0, 10000)
plt.ylim(0, 10000)
plt.title("Sample1")

a2 = fig.add_subplot(2, 2, 2)
a2.set_xscale('log')
a2.set_yscale('log')
plt.scatter(control[:-5], large_list[1][:-5],color='red')
plt.plot(range(1,10000),range(1,10000), 'k-')
plt.xlim(0, 10000)
plt.ylim(0, 10000)
plt.title("Sample2")

a3 = fig.add_subplot(2, 2, 3)
a3.set_xscale('log')
a3.set_yscale('log')
plt.scatter(control[:-5], large_list[2][:-5],color='orange')
plt.plot(range(1,10000),range(1,10000), 'k-')
plt.xlim(0, 10000)
plt.ylim(0, 10000)
plt.title("Sample3")

a4 = fig.add_subplot(2, 2, 4)
a4.set_xscale('log')
a4.set_yscale('log')
plt.scatter(control[:-5], large_list[3][:-5],color='green')
plt.plot(range(1,10000),range(1,10000), 'k-')
plt.xlim(0, 10000)
plt.ylim(0, 10000)
plt.title("Sample4")

plt.show()

sys.exit() 

out_list = []
count = 0
for gene in genes_of_interest :
    small_list = []
    small_list.append(gene)
    for a0line in a0lines :
        if a0line.split("\t")[0] == gene :
            small_list.append(a0line.split("\t")[1].strip("\n"))
    for large in large_list :
        for lar in large :
            if lar.split("\t")[0] == gene :
                small_list.append(lar.split("\t")[1].strip("\n"))
    print count, small_list
    count += 1

# check length of all small_lists