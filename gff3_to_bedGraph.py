import os,sys
import datetime 

## save regions 
type = sys.argv[1]
outfile1 = open(str(type)+"_BSseqcap_B73_Mo17_by_region.out", 'w')
outfile2 = open(str(type)+"_BSseqcap_B73_Mo17_by_cytosine.out", 'w')
outfile1.write("chr\tstart\tend\tB73 metC\tB73 unmetC\tMo17 metC\tMo17 unmetC\tOh43 metC\tOh43 unmetC\tTx303 metC\tTx303 unmetC\tCML322 metC\tCML322 unmetC\n")
outfile2.write("chr\tstart\tend\tB73 metC\tB73 unmetC\tMo17 metC\tMo17 unmetC\tOh43 metC\tOh43 unmetC\tTx303 metC\tTx303 unmetC\tCML322 metC\tCML322 unmetC\n")
for chr in range(1,11) :
    print datetime.datetime.now(), chr
    large_list = []
    for x in range(1,301354115):
        large_list.append([0]*10)
    print datetime.datetime.now(), "B73"
    for b in open("B73_all3/"+str(type)+"_context_B73_all3_R1_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(b.split("\t")[1])][0] += int(b.split("\t")[4])
        large_list[int(b.split("\t")[1])][1] += int(b.split("\t")[5].strip("\n"))
    print datetime.datetime.now(), "Mo17"
    for c in open("Mo17_all3/"+str(type)+"_context_Mo17_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(c.split("\t")[1])][2] += int(c.split("\t")[4])
        large_list[int(c.split("\t")[1])][3] += int(c.split("\t")[5].strip("\n"))
    print datetime.datetime.now(), "Oh43"
    for d in open("Oh43_all3/"+str(type)+"_context_Oh43_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(d.split("\t")[1])][4] += int(d.split("\t")[4])
        large_list[int(d.split("\t")[1])][5] += int(d.split("\t")[5].strip("\n"))
    print datetime.datetime.now(), "Tx303"
    for e in open("Tx303_all3/"+str(type)+"_context_Tx303_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(e.split("\t")[1])][6] += int(e.split("\t")[4])
        large_list[int(e.split("\t")[1])][7] += int(e.split("\t")[5].strip("\n"))
    print datetime.datetime.now(), "CML322"
    for f in open("CML322_all3/"+str(type)+"_context_CML322_all3_bt202_sortedn_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(f.split("\t")[1])][8] += int(f.split("\t")[4])
        large_list[int(f.split("\t")[1])][9] += int(f.split("\t")[5].strip("\n"))
    print datetime.datetime.now()
    regions = []
    for a in open("../BSseqcap_V1.sorted.gff3", 'r') :
        if a[0] != "#" and a.split("\t")[0][3:] == str(chr) :
            sum_met_1 = 0
            sum_unmet_1 = 0
            sum_met_2 = 0
            sum_unmet_2 = 0
            sum_met_3 = 0
            sum_unmet_3 = 0
            sum_met_4 = 0
            sum_unmet_4 = 0
            sum_met_5 = 0
            sum_unmet_5 = 0
            count = 0
            for vals in large_list[int(a.split("\t")[3]):int(a.split("\t")[4])+1] :
                sum_met_1 += vals[0]
                sum_unmet_1 += vals[1]
                sum_met_2 += vals[2]
                sum_unmet_2 += vals[3]
                sum_met_3 += vals[4]
                sum_unmet_3 += vals[5]
                sum_met_4 += vals[6]
                sum_unmet_4 += vals[7]
                sum_met_5 += vals[8]
                sum_unmet_5 += vals[9]
                outfile2.write("chr"+str(chr)+"\t"+str(int(a.split("\t")[3])+count)+"\t"+str(int(a.split("\t")[3])+count)+"\t"+str(vals[0])+"\t"+str(vals[1])+"\t"+str(vals[2])+"\t"+str(vals[3])+"\t"+str(vals[4])+"\t"+str(vals[5])+"\t"+str(vals[6])+"\t"+str(vals[7])+"\t"+str(vals[8])+"\t"+str(vals[9])+"\n")
                count += 1
            outfile1.write(a.split("\t")[0]+"\t"+a.split("\t")[3]+"\t"+a.split("\t")[4]+"\t"+str(sum_met_1)+"\t"+str(sum_unmet_1)+"\t"+str(sum_met_2)+"\t"+str(sum_unmet_2)+"\t"+str(sum_met_3)+"\t"+str(sum_unmet_3)+"\t"+str(sum_met_4)+"\t"+str(sum_unmet_4)+"\t"+str(sum_met_5)+"\t"+str(sum_unmet_5)+"\n")
outfile1.close()
outfile2.close()