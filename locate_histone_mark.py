### distribution of H3K4me3 mapping within gene

import os,sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

### make gene bed from gff3
outfile = open("Zea_mays.AGPv3.23.gene.gff3", 'w')
for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gff3", 'r') :
    if a[0] != "#" and a.split("\t")[2] == "gene" :
        outfile.write(a)
outfile.close()

peak_percentage_list = []
for chr in range(1,11) :
    print "Running chr"+str(chr)
    peaks_list = []
    for a in open(sys.argv[1], 'r') :
        if a.split("\t")[0] == str(chr) :
            peaks_list.append(int(float(int(a.split("\t")[1])+int(a.split("\t")[2]))/2))
    for a in open("Zea_mays.AGPv3.23.gene.gff3", 'r') :
        if a.split("\t")[0] == str(chr) :
            if a.split("\t")[6] == "+" :
                for peak in peaks_list :
                    ### within peak : calculate percentile
                    if int(peak) >= int(a.split("\t")[3]) and int(peak) <= int(a.split("\t")[4]) :
                        gene_length = float(int(a.split("\t")[4])-int(a.split("\t")[3]))
                        peak_location = float(peak-int(a.split("\t")[3]))
                        peak_percentage_list.append(peak_location/gene_length*100)
                    elif int(peak) >= int(a.split("\t")[3])-1000 and int(peak) < int(a.split("\t")[3]) :
                        peak_location = float(int(a.split("\t")[3])-peak)
                        peak_percentage_list.append((peak_location/10)*-1)
                    elif int(peak) > int(a.split("\t")[4]) and int(peak) <= int(a.split("\t")[4])+1000 :
                        peak_location = peak-int(a.split("\t")[4])
                        peak_percentage_list.append((peak_location/10)+100)
                    
            elif a.split("\t")[6] == "-" :
                for peak in peaks_list :
                    ### within peak : calculate percentile
                    if int(peak) >= int(a.split("\t")[3]) and int(peak) <= int(a.split("\t")[4]) :
                        gene_length = float(int(a.split("\t")[4])-int(a.split("\t")[3]))
                        peak_location = float(peak-int(a.split("\t")[3]))
                        peak_percentage_list.append(100-(peak_location/gene_length*100))
                    elif int(peak) <= int(a.split("\t")[4])+1000 and int(peak) > int(a.split("\t")[4]) :
                        peak_location = float(peak-int(a.split("\t")[4]))
                        peak_percentage_list.append((peak_location/10)*-1)
                    elif int(peak) < int(a.split("\t")[3]) and int(peak) >= int(a.split("\t")[3])-1000 :
                        peak_location = float(int(a.split("\t")[3])-peak)
                        peak_percentage_list.append((peak_location/10)+100)              
            
print len(peak_percentage_list)

plt.hist(peak_percentage_list, 50, facecolor='green', alpha=0.75)
plt.xlabel('peak location')
plt.ylabel('frequency')
plt.savefig('chip_location_'+str(sys.argv[1].split("."))+'.png')
                             