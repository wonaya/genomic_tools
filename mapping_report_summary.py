## for all bismark output returns mapping efficiency in percentage 
## by Jawon Song 01 March 2013
## Last edited 05 June 2013

import os,sys

outfile = open("mapping_efficiency_summary.txt", 'w')
for files in os.listdir("."):
    if files.endswith("bt2_paired-end_mapping_report.txt"):
        for line in open(files, 'r') :
            if line.split(":")[0] == "Mapping efficiency" :
                #outfile.write(str(files.split(".fastq")[0])+"\t"+str(line.split(":\t")[1]))
                outfile.write(str(files.split(".fastq")[0])+"\t"+str(line.split(":\t")[1]))
outfile.close()

os.system("sort -nk1 mapping_efficiency_summary.txt > mapping_efficiency_summary_sort.txt")
os.system("rm -Rf mapping_efficiency_summary.txt")