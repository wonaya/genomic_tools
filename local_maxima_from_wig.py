import os,sys
import numpy as np
import random
### find local maxima and take some window size around it to find enrichment vs. permuted region (same length) 
### promoters (upstream 500bp of TSS) to see enrichment vs. permuted region
"""
### save per chromosome in array
outfile  = open("ES_logFC.2.smooth.2C_broad_peaks.overlap.bedGraph", 'w')
#outfile2 = open("Random.1kb.2C_broad_peaks.2.overlap.bedGraph", 'w')
for chr in range(10,11) :
    print "chr"+str(chr)
    val_list = []
    loc_list = []
    for a in open("/scratch/02114/wonaya/NCSU_HiSeq/7-21-14_RelicationTiming_Hiseq/merged_1kb/ES_logFC.2.smooth.wig", 'r') :
        if a.split("\t")[0] == str(chr) : 
            val_list.append(float(a.split("\t")[3].strip("\n")))
            loc_list.append(int(a.split("\t")[1]))
    a = np.array(val_list, dtype=np.float)
    b = np.array(loc_list, dtype=np.int)
    minima_list = []
    count = 0
    for val in a :
        if float(val) > 0 : 
            if a[count+1] <= a[count] and a[count-1] <= a[count] :
                minima_list.append([a[count], loc_list[count]])
        count += 1
    
    ## read peak data
    peak_list = []
    for b in open("/scratch/02114/wonaya/NCSU_HiSeq/06-24-14_H3K56Ac_Hiseq/2C_broad_peaks.bed", 'r') :
        ### MACS score set as 100
        if b.split("\t")[0] == str(chr) and len(b.split("\t")) > 7 and int(b.split("\t")[4]) > 100 :
            peak_list.append([int(b.split("\t")[1]),int(b.split("\t")[2]),int(b.split("\t")[4])])
    
    for minima in minima_list :
        for peak in peak_list :
            ## window size set as 1kb
            ## if theres even 1bp overlap, call it
            if len(set(range(minima[1],minima[1]+1000))&set(range(peak[0],peak[1]))) > 0 :
                outfile.write(str(chr)+"\t"+str(minima[1])+"\t"+str(minima[1]+1000)+"\t"+str(peak[2])+"\n")
    
    round = 1
    while round <= 100 :
        outfile2 = open("Random.1kb.2C_broad_peaks."+str(round)+".overlap.bedGraph", 'w')
        print "permutation round:", round
        random_list = []
        x = 0
        while x < 881 : 
            value = random.choice(loc_list)
            for peak in peak_list :
                if len(set(range(int(value), int(value)+1000))&set(range(peak[0], peak[1]))) > 0 :
                    outfile2.write(str(chr)+"\t"+str(value)+"\t"+str(value+1000)+"\t"+str(peak[2])+"\n")
                    x += 1
        round += 1
        outfile2.close()
outfile.close() 
"""
"""
### question: how many times is peak found associated with local maxima
### question: is K56 peaks associated with Early replicating regions?

### find region with +ve replication
time = sys.argv[1]

outfile = open("MitS_"+str(time)+"_4x_sig_region.2.bedGraph",'w')
for chr in range(1,11) :
    val_list = []
    print "chr"+str(chr)
    for a in open(str(time)+".4x_logFC.2.smooth.wig", 'r') :
    #for a in open("Endo-"+str(time)+"_logFC.2.smooth.wig", 'r') :
        if a.split("\t")[0] == str(chr) : 
            val_list.append(float(a.split("\t")[3].strip("\n")))
    a = np.array(val_list, dtype=np.float)
    index = 0
    loc_list = []
    sig_range_list = []
    for val in a : 
        if index >= 1 :
            if val <= 0 :
                if a[index-1] > 0 :
                   loc_list.append(index)
                   if len(loc_list) > 1 :
                       if max(val_list[loc_list[0]:loc_list[1]]) >= 1 :
                           sig_range_list.append(loc_list)
                   loc_list = []
            elif val > 0 :
                if a[index-1] <= 0 :
                   loc_list.append(index)
        index += 1
    for sig_loc in sig_range_list :
        outfile.write(str(chr)+"\t")
        outfile.write(str(sig_loc[0]*1000)+"\t")
        outfile.write(str(sig_loc[1]*1000-1)+"\t")
        outfile.write("1\n")
outfile.close()

### draw broad peaks score distribution
"""
"""
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

x = []
for a in open("2C_broad_peaks.bed", 'r') :
    if len(a.split("\t")) > 4 : 
        x.append(float(a.split("\t")[4]))
print len(x)

numBins = 50
ax.hist(x,numBins,color='green',alpha=0.8)
plt.show()
plt.savefig("test.png")
"""
"""
### broad peaks bed to igv readable bed
time = sys.argv[1]

outfile = open(time+"_broad_peaks_S50.bedGraph", 'w') 
for a in open(time+"_broad_peaks.bed", 'r') :
    if int(a.split("\t")[4]) >= 50 : 
        outfile.write("\t".join(a.split("\t")[:3]))
        outfile.write("\t")
        outfile.write(a.split("\t")[4])
        outfile.write("\n")
outfile.close()
"""
"""
### open K56 data
### use intersectBed

### get size of each chromosome
genome_list = []
whole_list = []
chr_list = []
print "reading in genome"
for genome in open("/work/02114/wonaya/genome/Zea_mays.AGPv2.14.fa", 'r') :
    if genome[0] == ">" :
        if len(chr_list) == 0 :
            chr_list.append(genome.split(" dna")[0][1:])
        else :
            whole_list.append("".join(genome_list))
            chr_list.append(genome.split(" dna")[0][1:])
            genome_list = []
    else :
        genome_list.append(genome.strip("\n"))
whole_list.append("".join(genome_list))

length = 0
for chr in range(1,11) :
    length += len(whole_list[chr-1])
print "total Genome length", length

length = 0
for chr in range(1,11) :
    print chr, len(whole_list[chr-1])
    for a in open("ES_sig_region.bedGraph",'r') :
        if a.split("\t")[0] == str(chr) :
            length += int(a.split("\t")[2])-int(a.split("\t")[1])
print "ES", length

length = 0
for chr in range(1,11) :
    print chr, len(whole_list[chr-1])
    for a in open("MS_sig_region.bedGraph",'r') :
        if a.split("\t")[0] == str(chr) :
            length += int(a.split("\t")[2])-int(a.split("\t")[1])
print "MS", length

length = 0
for chr in range(1,11) :
    print chr, len(whole_list[chr-1])
    for a in open("LS_sig_region.bedGraph",'r') :
        if a.split("\t")[0] == str(chr) :
            length += int(a.split("\t")[2])-int(a.split("\t")[1])
print "LS", length
"""

### open RNAseq data
### draw distribution of FPKM data
"""
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

x = []
genes = []
for a in open("cufflinks_out/genes.fpkm_tracking", 'r') :
    if len(a.split("\t")) > 9 and a.split("\t")[0] != 'tracking_id' and float(a.split("\t")[9]) >= 10000 : 
        x.append(float(a.split("\t")[9]))
        genes.append(a.split("\t")[0])
print len(x)
print genes
numBins = 50
ax.hist(x,numBins,color='green',alpha=0.8)
plt.yscale('log',nonposy='clip')
plt.show()
plt.savefig("test.png")
"""
"""
### get mapman annotation of highest expressed genes
x = []
genes = []
for a in open("cufflinks_out/genes.fpkm_tracking", 'r') :
    if len(a.split("\t")) > 9 and a.split("\t")[0] != 'tracking_id' and float(a.split("\t")[9]) >= 1000 : 
        x.append(float(a.split("\t")[9]))
        genes.append(a.split("\t")[0])

for gene in genes :
    for a in open("/work/02114/wonaya/genome/annotation/MapMan_B73_annot.txt", 'r') :
        if a.split("\t")[0] == gene : 
            print gene, a
            
"""          
time = sys.argv[1]   
### get size 
total = 0
for a in open("MitS_"+str(time)+"_4x_sig_region.2.bedGraph", 'r') :
    total += int(a.split("\t")[2])-int(a.split("\t")[1])
print total
