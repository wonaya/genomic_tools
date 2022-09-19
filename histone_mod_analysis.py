import os,sys

## get distribution of p-values

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

x = []
for a in open(sys.argv[1], 'r') :
    x.append(float(a.split("\t")[7]))
plt.hist(x, 50, facecolor='green', alpha=0.75)
plt.show()

## get rid of all under -log(Pvalue) of 50

outfile = open(sys.argv[1].split(".broadPeak")[0]+"_S50.bedGraph", 'w')
for a in open(sys.argv[1], 'r') :
    if float(a.split("\t")[7]) >= 50 :
        outfile.write("\t".join(a.split("\t")[:3]))
        outfile.write("\t")
        outfile.write(a.split("\t")[7])
        outfile.write("\n")
outfile.close()

os.system("intersectBed -a ../06-24-14_H3K56Ac_Hiseq/ES_sig_region.bedGraph -b "+sys.argv[1].split(".broadPeak")[0]+"_S50.bedGraph > ES_"+sys.argv[1].split(".broadPeak")[0]+"_S50.bedGraph")
os.system("intersectBed -a ../06-24-14_H3K56Ac_Hiseq/MS_sig_region.bedGraph -b "+sys.argv[1].split(".broadPeak")[0]+"_S50.bedGraph > MS_"+sys.argv[1].split(".broadPeak")[0]+"_S50.bedGraph")
os.system("intersectBed -a ../06-24-14_H3K56Ac_Hiseq/LS_sig_region.bedGraph -b "+sys.argv[1].split(".broadPeak")[0]+"_S50.bedGraph > LS_"+sys.argv[1].split(".broadPeak")[0]+"_S50.bedGraph")
