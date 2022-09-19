import os,sys

for x in [408788,408790,408791,408794,408795,408796,408797,409086] :
    print "SRR"+str(x)
    os.system("samtools view -bS SRR"+str(x)+".fastq_bismark_bt2.sam > SRR"+str(x)+".fastq_bismark_bt2.bam")
    os.system("samtools flagstat SRR"+str(x)+".fastq_bismark_bt2.bam")