#pick 1000 regions at random 50, 100, 500
#calculate variability within all genos

import os,sys
import random
"""
for chr in range(2,11) :
    print chr
    outfile = open("HapMap_chr"+str(chr)+".bed", 'w')
    for a in open("/corral-repl/tacc/bio/jawon/HapMap_v3/HapMap_chr"+str(chr)+".txt", 'r') :
        if a.split("\t")[0] != "chrom" :
            outfile.write("\t".join(a.split("\t")[:2]))
            outfile.write("\t")
            outfile.write(str(float(len(a.split("\t")[2:])-(a.split("\t")[2:].count(str(a.split("\t")[2])+"\n")+a.split("\t")[2:].count(str(a.split("\t")[2]))))/float(len(a.split("\t")[2:]))*100))
            outfile.write("\n")
    outfile.close()
"""
print "reading Genome"
chr_list = []
count_list = []
count = 0
for genome in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.18.dna.fa", 'r') :
    if genome[0] == ">" :
        if count > 0 :
            count_list.append(count)
        chr_list.append(genome[1:].split(" ")[0].strip("\n"))
        count = 0
    else :
        count += len(genome)
count_list.append(count)
chr_dict = {}
for chr in chr_list :
    chr_dict[chr] = count_list[chr_list.index(chr)]


large_coordinate = []
large_hapval = []
large_cov = []
large_metrate = []
for chr in range(1,11) :
    print chr, "reading hapmap, metrate"
    coordinate = []
    hapval = []
    for hap in open("HapMap_chr"+str(chr)+".bed", 'r') :
        coordinate.append(int(hap.split("\t")[1]))
        hapval.append(float(hap.split("\t")[2].strip("\n")))
    large_coordinate.append(coordinate)
    large_hapval.append(hapval)

    cov = []
    metrate = []
    for covline in open("../10-01-13_5genos/B73_all3/CHH_context_B73_all3_R1_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        cov.append(int(covline.split("\t")[1]))
        metrate.append(float(covline.split("\t")[3]))
    large_cov.append(cov)
    large_metrate.append(metrate)
        
print "random sampling"    
outfile = open("snprate_vs_metrate.txt", 'w')
length = 500
x = 0
while x < 101 : 
    chr = random.randrange(1,11)
    max = chr_dict[str(chr)]
    coord = random.choice(large_coordinate[chr-1])
    if coord in large_cov[chr-1] :
        x += 1
        print x, chr, coord, large_hapval[chr-1][large_coordinate[chr-1].index(coord)], large_metrate[chr-1][large_cov[chr-1].index(coord)]
        outfile.write(str(chr)+"\t"+str(coord)+"\t"+str(large_hapval[chr-1][large_coordinate[chr-1].index(coord)])+"\t"+str(large_metrate[chr-1][large_cov[chr-1].index(coord)])+"\n")
outfile.close()
    