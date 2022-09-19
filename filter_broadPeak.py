import os,sys
import numpy as np
import matplotlib.pyplot as plt

### filter macs.broadPeak by percentile

score = []
for a in open(sys.argv[1], 'r') : 
    score.append(float(a.split("\t")[7]))

percentile_cutoff = float(np.percentile(score, 50))

print percentile_cutoff, len(score)

outfile = open(sys.argv[1].split(".")[0]+"_P50.broadPeak", 'w')
for a in open(sys.argv[1], 'r') :
    if float(a.split("\t")[7]) >= percentile_cutoff :
        outfile.write(a)
outfile.close()

plt.hist(score, 100, alpha=0.5, label='Peak', color='red')
plt.show() 