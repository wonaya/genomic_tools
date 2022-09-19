### ChIP data structure
import os,sys
import datetime

### EDIT input files
list_of_bams = [sys.argv[1],  sys.argv[2],sys.argv[3]]
input_bam = [sys.argv[4]]
fai = "/work/02114/wonaya/genome/Zea_mays.AGPv3.23/Zea_mays.AGPv3.23_without_scaffold.fa.fai"

outfile = open(sys.argv[5], 'w')

### EDIT comment section 
outfile.write("#author:jawon\n")
outfile.write("#date:"+str(datetime.date.today())+"\n")
outfile.write("#type:"+sys.argv[6]+"\n")
outfile.write("#comment:count_calculated_using_multiBamCov,Mapped_on_AGPv3.23,MACS2.1.0_used_to_get_peaks,3bioreps,50percent_overlap\n")
outfile.write("#columns:chr\tstart\tend\tsummed test count\tinput count\tMACS peak\tMACS score\n")

#### make 1kb window file
print "running make windows", datetime.datetime.now()
os.system("windowMaker -g "+str(fai)+" -w 1000 > bed_"+sys.argv[6]+".out")

### multibamcov
print "running multibamcov", datetime.datetime.now()
for a in input_bam+list_of_bams :
    if not os.path.isfile(a+".bai") :
        os.system("/opt/apps/samtools/0.1.19/samtools index "+a)
os.system("multiBamCov -bams "+" ".join(input_bam+list_of_bams)+" -bed bed_"+sys.argv[6]+".out > test_"+sys.argv[6]+".out")

### macs
print "running macs2", datetime.datetime.now()
os.system("macs2 callpeak --treatment "+" ".join(list_of_bams)+" --control "+" ".join(input_bam)+" -f BAM -g 2.3e9 --nomodel --nolambda --broad -n test_"+sys.argv[6])

### do parallel processing per chr
print "collecting data", datetime.datetime.now()
os.system("intersectBed -a test_"+sys.argv[6]+".out -b test_"+sys.argv[6]+"_peaks.broadPeak -f 0.50 -wa -wb > overlap_"+sys.argv[6]+".txt")

chr_dict = {}
diff_dict = {}
for a in open(str(fai), 'r') :
    chr_dict[a.split("\t")[0]] = []
    diff_dict[a.split("\t")[0]] = []
for a in open("overlap_"+sys.argv[6]+".txt", 'r') :
    chr_dict[a.split("\t")[0]].append(a.split("\t")[1]+"-"+a.split("\t")[2])
    diff_dict[a.split("\t")[0]].append(a.split("\t")[-2])

for a in open("test_"+sys.argv[6]+".out", 'r') :
    outfile.write("\t".join(a.split("\t")[:3])+"\t")
    total_test = sum(map(int, a.split("\t")[4:]))
    outfile.write(str(total_test)+"\t")
    outfile.write(a.split("\t")[3]+"\t")
    if a.split("\t")[1]+"-"+a.split("\t")[2] in chr_dict[a.split("\t")[0]] :
        index = chr_dict[a.split("\t")[0]].index(a.split("\t")[1]+"-"+a.split("\t")[2])
        score = diff_dict[a.split("\t")[0]][index]
        outfile.write("Yes\t"+str(score)+"\n")
    else :
        outfile.write("No\tNone\n")
outfile.close()
print "finishing up", datetime.datetime.now()

    
    