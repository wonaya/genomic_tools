### RNAseq data structure
import os,sys
import datetime

### chr start end count gene  cufflinks_fpkm
bam = [sys.argv[1]]
gff = "/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gff3"
fpkm = sys.argv[2]
fai = "/work/02114/wonaya/genome/Zea_mays.AGPv3.23/Zea_mays.AGPv3.23_without_scaffold.fa.fai"

outfile = open(sys.argv[3], 'w')

### EDIT comment section 
outfile.write("#author:jawon\n")
outfile.write("#date:"+str(datetime.date.today())+"\n")
outfile.write("#type:"+sys.argv[4]+"\n")
outfile.write("#comment:count_calculated_using_multiBamCov,mapped_on_AGPv3.23,Tophat_mapped,50percent_overlap\n")
outfile.write("#columns:chr\tstart\tend\tcount\tgene\tcufflinks_fpkm\n")

#### make 1kb window file
print "running make windows", datetime.datetime.now()
os.system("windowMaker -g "+str(fai)+" -w 1000 > bed_test.out")

### multibamcov
print "running multibamcov", datetime.datetime.now()
for a in bam :
    if not os.path.isfile(a+".bai") :
        os.system("/opt/apps/samtools/0.1.19/samtools index "+a)
os.system("multiBamCov -bams "+" ".join(bam)+" -bed bed_test.out > test_count.out")

### convert gff to bed
gff_file = open("gff_test.bed", 'w')
for a in open(gff, 'r') :
    if a[0] != "#" :
        if a.split("\t")[2] == "gene" :
            gff_file.write(a.split("\t")[0]+"\t"+a.split("\t")[3]+"\t"+a.split("\t")[4]+"\t")
            gene =  a.split("\t")[8].split(";")[0].split(":")[1]
            gff_file.write(gene+"\n")

### fpkm_dict
fpkm_dict = {}
fpkm_file = open(fpkm, 'r') 
fpkm_lines = fpkm_file.readlines()
for fpkm_line in fpkm_lines[1:]:
    fpkm_dict[fpkm_line.split("\t")[3]] = float(fpkm_line.split("\t")[9])

print "collecting data", datetime.datetime.now()
os.system("intersectBed -a test_count.out -b gff_test.bed -f 0.50 -wa -wb > overlap_test.txt")

chr_dict = {}
diff_dict = {}
for a in open(fai, 'r') :
    chr_dict[a.split("\t")[0]] = []
    diff_dict[a.split("\t")[0]] = []

for a in open("overlap_test.txt", 'r') :
    chr_dict[a.split("\t")[0]].append(a.split("\t")[1]+"-"+a.split("\t")[2])
    diff_dict[a.split("\t")[0]].append(a.split("\t")[-1].strip("\n"))

for a in open("test_count.out", 'r') :
    outfile.write("\t".join(a.split("\t")[:3])+"\t")
    outfile.write(a.split("\t")[3].strip("\n")+"\t")
    if a.split("\t")[1]+"-"+a.split("\t")[2] in chr_dict[a.split("\t")[0]] :
        index = chr_dict[a.split("\t")[0]].index(a.split("\t")[1]+"-"+a.split("\t")[2])
        score = diff_dict[a.split("\t")[0]][index]
        outfile.write("Yes\t"+str(score)+"\n")
    else :
        outfile.write("No\tNone\n")
outfile.close()
print "finishing up", datetime.datetime.now()

