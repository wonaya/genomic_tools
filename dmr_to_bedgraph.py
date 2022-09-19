import os,sys

outfile = open("CHG_DMR_5genos_merged_new2.bedGraph", 'w')
for a in open("CHG_DMR_5genos_merged_new2.txt", 'r') :
    outfile.write("\t".join(a.split("\t")[:3]))
    outfile.write("\t1\n")
outfile.close()

