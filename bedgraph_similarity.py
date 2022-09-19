#script.py first.bedgraph second.bedgraph

# sortBed -i Input_8C_allmerge_RPGCnorm_1kb.bedgraph > Input_8C_allmerge_RPGCnorm_1kb_sort.bedgraph
# intersectBed -wa -wb -a Input_8C_allmerge_RPGCnorm_1kb_sort.bedgraph -b maize_4.33_w1000.bed > Input_8C_allmerge_RPGCnorm_1kb_merge.bedgraph
# python merge_two_beds.py Input_2C_allmerge_RPGCnorm_1kb_merge.bedgraph Input_2C_allmerge_RPGCnorm_1kb_fix.bedgraph
# python similarity.py Input_2C_allmerge_RPGCnorm_1kb_fix.bedgraph Input_8C_allmerge_RPGCnorm_1kb_fix.bedgraph
import os,sys
import difflib
import numpy as np

os.system("/opt/apps/bedtools/2.25.0/bin/sortBed -i "+sys.argv[1]+" > "+sys.argv[1].split(".bedgraph")[0]+"_sort.bedgraph")
os.system("/opt/apps/bedtools/2.25.0/bin/sortBed -i "+sys.argv[2]+" > "+sys.argv[2].split(".bedgraph")[0]+"_sort.bedgraph")

sortfile1=sys.argv[1].split(".bedgraph")[0]+"_sort.bedgraph"
sortfile2=sys.argv[2].split(".bedgraph")[0]+"_sort.bedgraph"

os.system("/opt/apps/bedtools/2.25.0/bin/intersectBed -wa -wb -a "+sortfile1+" -b /work/02114/wonaya/genome/maize/Zea_mays.AGPv4.33/maize_4.33_w1000.bed > "+sortfile1.split("_sort.bedgraph")[0]+"_merge.bedgraph")
os.system("/opt/apps/bedtools/2.25.0/bin/intersectBed -wa -wb -a "+sortfile2+" -b /work/02114/wonaya/genome/maize/Zea_mays.AGPv4.33/maize_4.33_w1000.bed > "+sortfile2.split("_sort.bedgraph")[0]+"_merge.bedgraph")

mergefile1=sortfile1.split("_sort.bedgraph")[0]+"_merge.bedgraph"
mergefile2=sortfile2.split("_sort.bedgraph")[0]+"_merge.bedgraph"

file = open(mergefile1, 'r')
lines= file.readlines()
outfile1 = open(mergefile1.split("_merge.bedgraph")[0]+"_fix.bedgraph", 'w')
for a in lines :
    outfile1.write(a.split("\t")[-3])
    outfile1.write("\t")
    outfile1.write(a.split("\t")[-2])
    outfile1.write("\t")
    outfile1.write(a.split("\t")[-1].strip("\n"))
    outfile1.write("\t")
    outfile1.write(a.split("\t")[-4])
    outfile1.write("\n")
outfile1.close()

file = open(mergefile2, 'r')
lines= file.readlines()
outfile2 = open(mergefile2.split("_merge.bedgraph")[0]+"_fix.bedgraph", 'w')
for a in lines :
    outfile2.write(a.split("\t")[-3])
    outfile2.write("\t")
    outfile2.write(a.split("\t")[-2])
    outfile2.write("\t")
    outfile2.write(a.split("\t")[-1].strip("\n"))
    outfile2.write("\t")
    outfile2.write(a.split("\t")[-4])
    outfile2.write("\n")
outfile2.close()

fixfile1=mergefile1.split("_merge.bedgraph")[0]+"_fix.bedgraph"
fixfile2=mergefile2.split("_merge.bedgraph")[0]+"_fix.bedgraph"

list1 = []
for a in open(fixfile1, 'r'):
    list1.append(float(a.split("\t")[-1].strip("\n")))
list2 = []
for b in open(fixfile2, 'r'):
    list2.append(float(b.split("\t")[-1].strip("\n")))

print len(list1), len(list2)
#sm=difflib.SequenceMatcher(None, list1, list2)
#print sm.ratio()
print np.corrcoef(list1, list2)


