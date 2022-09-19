"""
### trim ###

import os,sys

outfile = open("trim_all.txt", 'w')
for ROOT,DIR,FILES in os.walk("."):
    for file in FILES:
       if file.endswith(('_R1_001.fastq')):
          file_id = "_".join(file.split("_")[:3])
          barcode = file.split("_")[1]
          print file_id, barcode
          outfile.write("/work/02114/wonaya/software/trim_galore --fastqc -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"+barcode+"ATCTCGTATGCCGTCTTCTGCTTG -a2 GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"+barcode+"ATCTCGTATGCCGTCTTCTGCTTG -t --paired "+file_id+"_R1_001.fastq "+file_id+"_R2_001.fastq")
          outfile.write("\n")
outfile.close()
"""
"""
### bwa mem ###

import os,sys

outfile = open("bwamem_all.txt", 'w')
for ROOT,DIR,FILES in os.walk("."):
    for file in FILES:
       if file.endswith(('_R1_001_val_1.fq')):
          file_id = "_".join(file.split("_")[:3])
          print file_id
          outfile.write("bwa mem /work/02114/wonaya/genome/Arabidopsis_thaliana.TAIR10.14/Arabidopsis_thaliana.TAIR10.14.dna.toplevel.fa "+file_id+"_R1_001_val_1.fq "+file_id+"_R2_001_val_2.fq > "+file_id+".sam")
          outfile.write("\n")
outfile.close()
"""
"""
### tophat ###
import os,sys

outfile = open("tophat_all.txt", 'w')
for ROOT,DIR,FILES in os.walk("."):
    for file in FILES:
       if file.endswith(('_R1_001_val_1.fq')):
          file_id = "_".join(file.split("_")[:3])
          print file_id
          outfile.write("tophat -o tophat_323_"+file_id+" -G /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gtf --no-novel-juncs /work/02114/wonaya/genome/Zea_mays.AGPv3.23/Zea_mays.AGPv3.23 "+file_id+"_R1_001_val_1.fq "+file_id+"_R2_001_val_2.fq")
          outfile.write("\n")
outfile.close()
"""
"""
### star ###
import os,sys

outfile = open("star_all.txt", 'w')
for ROOT,DIR,FILES in os.walk("."):
    for file in FILES:
       if file.endswith(('_R1_001_val_1.fq')):
          file_id = "_".join(file.split("_")[:3])
          print file_id
          outfile.write("/work/02114/wonaya/software/STAR-STAR_2.4.0f1/bin/Linux_x86_64/STAR --outSAMstrandField intronMotif --runThreadN 12 --genomeDir /work/02114/wonaya/genome/Zea_mays.AGPv3.23/ --readFilesIn "+file_id+"_R1_001_val_1.fq "+file_id+"_R2_001_val_2.fq --outFileNamePrefix STAR_"+file_id)
          outfile.write("\n")
outfile.close()
"""
"""
### sam to sorted bam ###
import os,sys

outfile = open("sam2bam.txt", 'w')
for ROOT,DIR,FILES in os.walk("."):
    for file in FILES:
       if file.endswith(('.sam')):
          file_id = file.split(".")[0]
          print file_id
          outfile.write("samtools view -bS "+file_id+".sam | samtools sort - "+file_id+".sorted")
          outfile.write("\n")
outfile.close()
"""
### cufflinks ###
import os,sys
outfile = open("cufflinks.txt", 'w')
for ROOT,DIR,FILES in os.walk("."):
    for file in FILES:
       if file.endswith(('L002Aligned.sorted.bam')):
          file_id = file.split(".")[0]
          print file_id
          outfile.write("cufflinks -o "+file_id+" -G /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gtf "+file)
          outfile.write("\n")
outfile.close()
  
"""
### picard rmdup ###
import os,sys

outfile = open("picard_rmdup.txt", 'w')
for ROOT,DIR,FILES in os.walk("."):
    for file in FILES:
       if file.endswith(('.bam.sorted.bam')):
           file_id = file.split(".bam")[0]
           print file_id
           outfile.write("java -Xmx4g -jar /opt/apps/picard/1.107/MarkDuplicates.jar INPUT="+file_id+".bam.sorted.bam OUTPUT="+file_id+".rmdup.bam METRICS_FILE=dup.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true TMP_DIR=/tmp ASSUME_SORTED=true")
           outfile.write("\n")
outfile.close()
"""
"""
### picard downsample ###
import os,sys

outfile = open("picard_down.txt", 'w')
for ROOT,DIR,FILES in os.walk("."):
    for file in FILES:
       if file.endswith(('.rmdup.bam')):
           file_id = file.split(".rmdup")[0]
           print file_id
           outfile.write("java -Xmx4g -jar /opt/apps/picard/1.107/DownsampleSam.jar INPUT="+file_id+".rmdup.bam OUTPUT="+file_id+".001.ds.bam VALIDATION_STRINGENCY=SILENT P=0.01")
           outfile.write("\n")
           outfile.write("java -Xmx4g -jar /opt/apps/picard/1.107/DownsampleSam.jar INPUT="+file_id+".rmdup.bam OUTPUT="+file_id+".01.ds.bam VALIDATION_STRINGENCY=SILENT P=0.1")
           outfile.write("\n")
           outfile.write("java -Xmx4g -jar /opt/apps/picard/1.107/DownsampleSam.jar INPUT="+file_id+".rmdup.bam OUTPUT="+file_id+".05.ds.bam VALIDATION_STRINGENCY=SILENT P=0.5")
           outfile.write("\n")
outfile.close()
"""