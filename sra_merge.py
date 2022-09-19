import os,sys

"""
## merge
for no in range(7,93) :
    id = str(521000+no)
    print id
    os.system("/opt/apps/samtools/0.1.19/samtools view -bS SRR"+str(id)+".fastq_bismark_bt2.sam > SRR"+str(id)+".fastq_bismark_bt2.bam")

list = []
for no in range(7,50) :
    list.append("SRR"+str(id)+".fastq_bismark_bt2.bam")
os.system("/opt/apps/samtools/0.1.19/samtools merge B73_07to49.bam "+" ".join(list))
os.system("/opt/apps/samtools/0.1.19/samtools view -h -o B73_07to49.sam B73_07to49.bam")

list = []
for no in range(50,93) :
    list.append("SRR"+str(id)+".fastq_bismark_bt2.bam")
os.system("/opt/apps/samtools/0.1.19/samtools merge Mo17_50to93.bam "+" ".join(list))
os.system("/opt/apps/samtools/0.1.19/samtools view -h -o Mo17_50to93.sam Mo17_50to93.bam")

os.system("/opt/apps/bismark/0.10.1/bismark_methylation_extractor -s --report --bedGraph B73_07to49.sam")
os.system("/opt/apps/bismark/0.10.1/bismark_methylation_extractor -s --report --bedGraph Mo17_50to93.sam")

for type in ['CpG','CHG','CHH'] :
    print "B73", type
    otfile = open(type+"_OT_B73.txt",'w')
    obfile = open(type+"_OB_B73.txt",'w')
    ctotfile = open(type+"_CTOT_B73.txt",'w')
    ctobfile = open(type+"_CTOB_B73.txt",'w')
    for no in range(7,50) :
        id = no+521000
        print id
        a1 = open(type+"_OT_SRR"+str(id)+".fastq_bismark_bt2.txt", 'r') 
        a1lines = a1.readlines()
        for a1line in a1lines[1:] :
            otfile.write(a1line)
        a2 = open(type+"_OB_SRR"+str(id)+".fastq_bismark_bt2.txt", 'r') 
        a2lines = a2.readlines()
        for a2line in a2lines[1:] :
            obfile.write(a2line)
        a3 = open(type+"_CTOT_SRR"+str(id)+".fastq_bismark_bt2.txt", 'r') 
        a3lines = a3.readlines()
        for a3line in a3lines[1:] :
            ctotfile.write(a3line)
        a4 = open(type+"_CTOB_SRR"+str(id)+".fastq_bismark_bt2.txt", 'r') 
        a4lines = a4.readlines()
        for a4line in a4lines[1:] :
            ctobfile.write(a4line)
    otfile.close()
    obfile.close()
    ctotfile.close()
    ctobfile.close()
"""
for type in ['CpG','CHG','CHH'] :
    print "B73", type, "merge"
    otfile = open(type+"_OT_B73_dawe.txt",'w')
    obfile = open(type+"_OB_B73_dawe.txt",'w')
    ctotfile = open(type+"_CTOT_B73_dawe.txt",'w')
    ctobfile = open(type+"_CTOB_B73_dawe.txt",'w')
    for id in [8788,8790,8791,8794,8795,8796,8797,9086] :
        print id
        a1 = open(type+"_OT_SRR40"+str(id)+".fastq_bismark_bt2.txt", 'r') 
        a1lines = a1.readlines()
        for a1line in a1lines[1:] :
            otfile.write(a1line)
        a2 = open(type+"_OB_SRR40"+str(id)+".fastq_bismark_bt2.txt", 'r') 
        a2lines = a2.readlines()
        for a2line in a2lines[1:] :
            obfile.write(a2line)
        a3 = open(type+"_CTOT_SRR40"+str(id)+".fastq_bismark_bt2.txt", 'r') 
        a3lines = a3.readlines()
        for a3line in a3lines[1:] :
            ctotfile.write(a3line)
        a4 = open(type+"_CTOB_SRR40"+str(id)+".fastq_bismark_bt2.txt", 'r') 
        a4lines = a4.readlines()
        for a4line in a4lines[1:] :
            ctobfile.write(a4line)
    otfile.close()
    obfile.close()
    ctotfile.close()
    ctobfile.close()
    context = type
    print "B73", context, "readin"
    large_list = []
    for x in range(1,150000000) :
        large_list.append([0]*4)

    for a in open(context+"_OT_B73_dawe.txt", 'r') :
        if a.split("\t")[2] == "10" :
            if a.split("\t")[1] == "+" :
                large_list[int(a.split("\t")[3])][0] += 1
            elif a.split("\t")[1] == "-" : 
                large_list[int(a.split("\t")[3])][1] += 1
    for a in open(context+"_OB_B73_dawe.txt", 'r') :
        if a.split("\t")[2] == "10" :
            if a.split("\t")[1] == "+" :
                large_list[int(a.split("\t")[3])][0] += 1
            elif a.split("\t")[1] == "-" : 
                large_list[int(a.split("\t")[3])][1] += 1
    print "B73", type, "bedGraph"
    index = 0
    outfile = open(context+"_B73_dawe_chr10.bedGraph",'w')
    for list in large_list :
        if list[:2] != [0]*2 :
            outfile.write("10\t")
            outfile.write(str(index-1)+"\t"+str(index)+"\t")
            outfile.write(str((float(list[0])/(float(list[0])+float(list[1]))*100)))
            outfile.write("\n")
        index += 1
    outfile.close()

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
    outfile.write("chr\tstart\tend\tB73 met%\tRob B73 met count\tRob B73 unmet count\tRob B73 met%\n")
    for aline in alines[1:]:
        ranges = range(int(aline.split("\t")[1]),int(aline.split("\t")[2])+1)
        if len(ranges) >= 50 : ## if DMR length >= 50
            b73_cg = (aline.split("\t")[fromto[0]])
            b73_rob_met = 0
            b73_rob_unmet = 0
            for count in large_list[int(aline.split("\t")[1]):int(aline.split("\t")[2])+1] :
                b73_rob_met += count[0]
                b73_rob_unmet += count[1]
            if float(b73_rob_met)+float(b73_rob_unmet) <= 0 : 
                b73_rob_metperc = "NA"
            else :
                b73_rob_metperc = float(b73_rob_met)/(float(b73_rob_met)+float(b73_rob_unmet))
            outfile.write(str("\t".join(aline.split("\t")[:3]))+"\t"+str(b73_cg)+"\t"+str(b73_rob_met)+"\t"+str(b73_rob_unmet)+"\t"+str(b73_rob_metperc)+"\n")
            del b73_rob_metperc
    outfile.close()        