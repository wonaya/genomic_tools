import os,sys
import datetime

print datetime.datetime.now()
large_list = []
for x in range(1,150000000) :
    large_list.append([0]*4)

context = sys.argv[1]

for a in open("GSE39232_B73_sequence.rmdup.chr10.bed", 'r') :
    if a.split("\t")[9].strip("\n") == context :
        large_list[int(a.split("\t")[2])][0] += int(a.split("\t")[6])
        large_list[int(a.split("\t")[2])][1] += int(a.split("\t")[7])
print "b73 done"
for a in open("GSE39232_MO17_sequence.rmdup.chr10.bed", 'r') :
    if a.split("\t")[9].strip("\n") == context :
        large_list[int(a.split("\t")[2])][2] += int(a.split("\t")[6])
        large_list[int(a.split("\t")[2])][3] += int(a.split("\t")[7])
print "mo17 done"

a = open("../10-01-13_5genos/"+str(context)+"_DMR_5genos_chr10_new2.txt", 'r')
dmr_coord = [] 
alines = a.readlines()

if context == "CpG" : 
    fromto = [10,11]
elif context == "CHG" :
    fromto = [15,16]
elif context == "CHH" :
    fromto = [20,21]

outfile = open(str(context)+"_DMR_metval_compare_chr10.txt", 'w')
outfile.write("chr\tstart\tend\tB73 met%\tMo17 met%\tRob B73 met count\tRob B73 unmet count\tRob Mo17 met count\tRob Mo17 unmet count\tRob B73 met%\tRob Mo17 met%\n")
for aline in alines[1:]:
    ranges = range(int(aline.split("\t")[1]),int(aline.split("\t")[2])+1)
    if len(ranges) >= 50 : ## if DMR length >= 50
        b73_cg = (aline.split("\t")[fromto[0]])
        mo17_cg = (aline.split("\t")[fromto[1]])
        b73_rob_met = 0
        b73_rob_unmet = 0
        mo17_rob_met = 0
        mo17_rob_unmet = 0
        for count in large_list[int(aline.split("\t")[1]):int(aline.split("\t")[2])+1] :
            b73_rob_met += count[0]
            b73_rob_unmet += count[1]
            mo17_rob_met += count[2]
            mo17_rob_unmet += count[3]
        if float(b73_rob_met)+float(b73_rob_unmet) <= 0 : 
            b73_rob_metperc = "NA"
        else :
            b73_rob_metperc = float(b73_rob_met)/(float(b73_rob_met)+float(b73_rob_unmet))
        if float(mo17_rob_met)+float(mo17_rob_unmet) <= 0 : 
            mo17_rob_metperc = "NA"
        else: 
            mo17_rob_metperc = float(mo17_rob_met)/(float(mo17_rob_met)+float(mo17_rob_unmet))   
        print str("\t".join(aline.split("\t")[:3]))+"\t"+str(b73_cg)+"\t"+str(mo17_cg)+"\t"+str(b73_rob_met)+"\t"+str(b73_rob_unmet)+"\t"+str(mo17_rob_met)+"\t"+str(mo17_rob_unmet)+"\t"+str(b73_rob_metperc)+"\t"+str(mo17_rob_metperc)+"\n"
        outfile.write(str("\t".join(aline.split("\t")[:3]))+"\t"+str(b73_cg)+"\t"+str(mo17_cg)+"\t"+str(b73_rob_met)+"\t"+str(b73_rob_unmet)+"\t"+str(mo17_rob_met)+"\t"+str(mo17_rob_unmet)+"\t"+str(b73_rob_metperc)+"\t"+str(mo17_rob_metperc)+"\n")
        del b73_rob_metperc
        del mo17_rob_metperc
outfile.close()
print datetime.datetime.now()
