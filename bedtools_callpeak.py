## module load samtools bedtools

import os,sys

#################################
#### PART 1 #####################
#################################

#### Sample 1
### convert test1 bam to bed using deeptools
os.system("samtools index test1.bam")
os.system("bamCoverage --bam test1.bam --outFileFormat bedgraph -o test1_bin1kb.bed --binSize 1000 >> tmp.txt")
os.system("sortBed -i test1_bin1kb.bed > test1_bin1kb_sort.bed")
os.system("intersectBed -wa -wb -a test1_bin1kb_sort.bed -b maize_4.33_w1000.bed > test1_bin1kb_fix.bedgraph")

## convert control1 bam to bed using deeptools
os.system("samtools index control1.bam")
os.system("bamCoverage --bam control1.bam --outFileFormat bedgraph -o control1_bin1kb.bed --binSize 1000 >> tmp.txt")
os.system("sortBed -i control1_bin1kb.bed > control1_bin1kb_sort.bed")
os.system("intersectBed -wa -wb -a control1_bin1kb_sort.bed -b maize_4.33_w1000.bed > control1_bin1kb_fix.bedgraph")

## merge two test1 bed files into 1
file = open('test1_bin1kb_fix.bedgraph', 'r')
lines= file.readlines()
outfile = open('test1_bin1kb_merge.bedgraph', 'w')
for a in lines :
    outfile.write(a.split("\t")[-3])
    outfile.write("\t")
    outfile.write(a.split("\t")[-2])
    outfile.write("\t")
    outfile.write(a.split("\t")[-1].strip("\n"))
    outfile.write("\t")
    outfile.write(a.split("\t")[-4])
    outfile.write("\n")
outfile.close()

file = open('control1_bin1kb_fix.bedgraph', 'r')
lines= file.readlines()
outfile = open('control1_bin1kb_merge.bedgraph', 'w')
for a in lines :
    outfile.write(a.split("\t")[-3])
    outfile.write("\t")
    outfile.write(a.split("\t")[-2])
    outfile.write("\t")
    outfile.write(a.split("\t")[-1].strip("\n"))
    outfile.write("\t")
    outfile.write(a.split("\t")[-4])
    outfile.write("\n")
outfile.close()

## at this point two files should be of exact length and first 3 columns should match exactly
## identify p<0.05 peaks using test1 bed 
val1_list = []
for a in open('test1_bin1kb_merge.bedgraph', 'r') :
    val1_list.append(float(a.split("\t")[-1].strip("\n")))
val2_list = []
for a in open('control1_bin1kb_merge.bedgraph', 'r') :
    val2_list.append(float(a.split("\t")[-1].strip("\n")))
diff_list = []
index = 0
for val1 in val1_list :
    diff_list.append(val1-val2_list[index])
    index += 1

print "test1", max(diff_list), min(diff_list) 

# histogram of this and save
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.hist(diff_list, 50, facecolor='green', alpha=0.75)
plt.savefig('hist_test1.png')
#plt.show()

#### Sample 2
### convert test2 bam to bed using deeptools
os.system("samtools index test2.bam")
os.system("bamCoverage --bam test2.bam --outFileFormat bedgraph -o test2_bin1kb.bed --binSize 1000 >> tmp.txt")
os.system("sortBed -i test2_bin1kb.bed > test2_bin1kb_sort.bed")
os.system("intersectBed -wa -wb -a test2_bin1kb_sort.bed -b maize_4.33_w1000.bed > test2_bin1kb_fix.bedgraph")

### convert control2 bam to bed using deeptools
os.system("samtools index control2.bam")
os.system("bamCoverage --bam control2.bam --outFileFormat bedgraph -o control2_bin1kb.bed --binSize 1000 >> tmp.txt")
os.system("sortBed -i control2_bin1kb.bed > control2_bin1kb_sort.bed")
os.system("intersectBed -wa -wb -a control2_bin1kb_sort.bed -b maize_4.33_w1000.bed > control2_bin1kb_fix.bedgraph")

### merge two test2 bed files into 1
file = open('test2_bin1kb_fix.bedgraph', 'r')
lines= file.readlines()
outfile = open('test2_bin1kb_merge.bedgraph', 'w')
for a in lines :
    outfile.write(a.split("\t")[-3])
    outfile.write("\t")
    outfile.write(a.split("\t")[-2])
    outfile.write("\t")
    outfile.write(a.split("\t")[-1].strip("\n"))
    outfile.write("\t")
    outfile.write(a.split("\t")[-4])
    outfile.write("\n")
outfile.close()

file = open('control2_bin1kb_fix.bedgraph', 'r')
lines= file.readlines()
outfile = open('control2_bin1kb_merge.bedgraph', 'w')
for a in lines :
    outfile.write(a.split("\t")[-3])
    outfile.write("\t")
    outfile.write(a.split("\t")[-2])
    outfile.write("\t")
    outfile.write(a.split("\t")[-1].strip("\n"))
    outfile.write("\t")
    outfile.write(a.split("\t")[-4])
    outfile.write("\n")
outfile.close()

## identify p<0.05 peaks using test2 bed 
val1_list = []
for a in open('test2_bin1kb_merge.bedgraph', 'r') :
    val1_list.append(float(a.split("\t")[-1].strip("\n")))
val2_list = []
for a in open('control2_bin1kb_merge.bedgraph', 'r') :
    val2_list.append(float(a.split("\t")[-1].strip("\n")))
diff_list = []
index = 0
for val1 in val1_list :
    diff_list.append(val1-val2_list[index])
    index += 1

print "test2", max(diff_list), min(diff_list) 

# histogram of this and save
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.hist(diff_list, 50, facecolor='green', alpha=0.75)
plt.savefig('hist_test2.png')
#plt.show()

os.system("rm -Rf *_fix.bedgraph *_bin1kb.bed *.bai")

#################################
#### PART 2 #####################
#################################

### identify p<0.05 peaks using bed 
low=-5000
high=5000

print "threshold", low, high
# identify regions that are 5% and 95% quartiles
file1 = open('test1_bin1kb_merge.bedgraph', 'r')
lines1 = file1.readlines()
file2 = open('control1_bin1kb_merge.bedgraph', 'r')
lines2 = file2.readlines()

outfile = open("diff_test1.bedgraph", 'w')
index = 0
for line1 in lines1 :
    diffval = float(line1.split("\t")[-1].strip("\n"))-float(lines2[index].split("\t")[-1].strip("\n"))
    if diffval < low or diffval > high :
        #print line1, lines2[index] ; sys.exit()
        outfile.write(line1.split("\t")[0]+"\t")
        outfile.write(line1.split("\t")[1]+"\t")
        outfile.write(line1.split("\t")[2]+"\t")
        outfile.write(str(diffval)+"\n")
outfile.close()

# identify regions that are 5% and 95% quartiles
file1 = open('test2_bin1kb_merge.bedgraph', 'r')
lines1 = file1.readlines()
file2 = open('control2_bin1kb_merge.bedgraph', 'r')
lines2 = file2.readlines()

outfile = open("diff_test2.bedgraph", 'w')
index = 0
for line1 in lines1 :
    diffval = float(line1.split("\t")[-1].strip("\n"))-float(lines2[index].split("\t")[-1].strip("\n"))
    if diffval < low or diffval > high :
        #print line1, lines2[index] ; sys.exit()
        outfile.write(line1.split("\t")[0]+"\t")
        outfile.write(line1.split("\t")[1]+"\t")
        outfile.write(line1.split("\t")[2]+"\t")
        outfile.write(str(diffval)+"\n")
outfile.close()

#################################
#### PART 3 #####################
#################################

