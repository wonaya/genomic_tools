### size distribution of segments

#!/usr/bin/env python
  
import os,sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
import datetime
from math import *
from scipy import stats
from optparse import OptionParser, OptionGroup
import pandas as pd
from urllib2 import urlopen
import math
import random

reptiming_gff3 = sys.argv[1]

### prepare dictionaries for reptiming
reptiming_dict = {}
reptiming_distance_dict = {}
for a in open(reptiming_gff3, 'r') :
    if a[0] != "#" :
        reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
        reptiming_distance_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]] = []
            
for a in open(reptiming_gff3, 'r') :
    if a[0] != "#" :
        chr = a.split("\t")[0]
        coord_start = int(a.split("\t")[3])-1
        coord_end   = int(a.split("\t")[4])-1
        reptiming_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]].append(chr+":"+str(coord_start)+"-"+str(coord_end))
        reptiming_distance_dict[a.split("\t")[8].split("Name=")[1].split(";")[0]].append(int(coord_end)-int(coord_start)+1)
### write into bed and sort
print np.percentile(reptiming_distance_dict['ES'], 5)
print np.percentile(reptiming_distance_dict['MS'], 5)
print np.percentile(reptiming_distance_dict['LS'], 5)
    
bins = np.linspace(0,100000, 100)

plt.hist(reptiming_distance_dict['ES'], bins, alpha=0.5, label='ES')
#plt.hist(reptiming_distance_dict['ESMS'], bins, alpha=0.5, label='ESMS')
plt.hist(reptiming_distance_dict['MS'], bins, alpha=0.5, label='MS')
#plt.hist(reptiming_distance_dict['MSLS'], bins, alpha=0.5, label='MSLS')
plt.hist(reptiming_distance_dict['LS'], bins, alpha=0.5, label='LS')
#plt.hist(reptiming_distance_dict['ESLS'], bins, alpha=0.5, label='ESLS')
#plt.hist(reptiming_distance_dict['ESMSLS'], bins, alpha=0.5, label='ESMSLS')
plt.legend(loc='upper right')
    
plt.savefig("distribution_distance.png")
    
        
        
