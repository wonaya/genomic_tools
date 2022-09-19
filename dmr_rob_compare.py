### 

# chr start end B73_metval Mo17_metval Rob_B73 Rob_Mo17 type
import os,sys
import datetime
#[1,2]: [met, unmet]
"""
outfile = open("GSE39232_B73_sequence.rmdup.chr10.bed", 'w')     
for a in open("GSE39232_B73_sequence.rmdup.bed", 'r') :
    if a.split("\t")[0] == "chr10" :
        outfile.write(a)
outfile.close()

print datetime.datetime.now()
outfile = open("GSE39232_MO17_sequence.rmdup.chr10.bed", 'w')     
for a in open("GSE39232_MO17_sequence.rmdup.bed", 'r') :
    if a.split("\t")[0] == "chr10" :
        outfile.write(a)
outfile.close()
"""
print "split file done", datetime.datetime.now()
a = open("../10-01-13_5genos/CpG_DMR_5genos_chr10_new2.txt", 'r')
dmr_coord = [] 
alines = a.readlines()
#print alines[0]
#print alines[0].split("\t")[10:12], alines[0].split("\t")[20:22] ### CHH methylation %
for aline in alines[1:]:
    ranges = range(int(aline.split("\t")[1]),int(aline.split("\t")[2])+1)
    if len(ranges) >= 50 : ## if DMR length >= 50
        b73_cg = (aline.split("\t")[10])
        mo17_cg = (aline.split("\t")[11])
        dmr_coord.append([ranges, b73_cg, mo17_cg,0,0,0,0])
        #b73_chh = (aline.split("\t")[20])
        #mo17_chh = (aline.split("\t")[21])
        #dmr_coord.append([ranges, b73_chh, mo17_chh,0,0,0,0])
print "making dmr list file done", datetime.datetime.now()

for a in open("GSE39232_B73_sequence.rmdup.chr10.bed", 'r') :
    for b in dmr_coord :
        if int(a.split("\t")[2]) in b[0] and a.split("\t")[9].strip("\n") == "CpG" :
        #if int(a.split("\t")[2]) in b[0] and a.split("\t")[9].strip("\n") == "CHH" :
            b[3] += int(a.split("\t")[6])
            b[4] += int(a.split("\t")[7])
print "making b73 summary done", datetime.datetime.now()
for a in open("GSE39232_MO17_sequence.rmdup.chr10.bed", 'r') :
    for b in dmr_coord :
        if int(a.split("\t")[2]) in b[0] and a.split("\t")[9].strip("\n") == "CpG" :
        #if int(a.split("\t")[2]) in b[0] and a.split("\t")[9].strip("\n") == "CHH":
            b[5] += int(a.split("\t")[6])
            b[6] += int(a.split("\t")[7])
print "making MO17 summary done",datetime.datetime.now()
 
outfile = open("Springer_DMR_GSE39232_CpG_comparison.txt", 'w')
outfile.write("chr\tstart\tend\t5geno B73 met\t5geno Mo17 met\tGSE B73 met\tGSE Mo17 met\n")           
for dmr in dmr_coord :
    outfile.write("chr10\t")
    outfile.write(dmr[0][0])
    outfile.write("\t")
    outfile.write(dmr[0][-1])
    outfile.write("\t")
    outfile.write(dmr[1])
    outfile.write("\t")
    outfile.write(dmr[2])
    outfile.write("\t")
    outfile.write(str(float(dmr[3])/(float(dmr[3])+float(dmr[4]))))
    outfile.write("\t")
    outfile.write(str(float(dmr[5])/(float(dmr[5])+float(dmr[6]))))
    outfile.write("\n")
outfile.close()
print "done"