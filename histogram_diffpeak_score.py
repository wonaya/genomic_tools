## histogram

import os,sys

## get distribution of p-values

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

x = []
for a in open(sys.argv[1], 'r') :
    if len(a.split("\t")) > 3 :
        x.append(float(a.split("\t")[4].strip("\n")))
plt.hist(x, 50, facecolor='green', alpha=0.75)
plt.show()
