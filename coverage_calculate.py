import os,sys
### coverage calculator

### usage coverage_calculate.py [bam] [fai] [-P] [read length]

### get genome size

genome_size = 0
for fai in open(sys.argv[2], 'r') :
    genome_size += int(fai.split("\t")[1])

os.system("samtools flagstat "+sys.argv[1]+" > flagstat.out")
flagstat = open("flagstat.out", 'r') 
a = flagstat.readlines()
mapped_reads = a[2].split(" + ")[0]
if sys.argv[3] == "-P" :
    total_base = int(mapped_reads)*2*int(sys.argv[4])
elif sys.argv[3] == "-S" :
    total_base = int(mapped_reads)*int(sys.argv[4])
coverage = float(total_base)/float(genome_size)

print coverage
    
    
