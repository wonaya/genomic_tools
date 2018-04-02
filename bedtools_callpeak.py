## module load samtools bedtools

import os,sys

### convert bam to bed using deeptools
samtools index test.bam
bamCoverage --bam test.bam --outFileFormat bedgraph -o test_bin10.bed --binSize 10
sortBed -i test_bin10.bed > test_bin10_sort.bed

samtools index control.bam 
bamCoverage --bam control.bam --outFileFormat bedgraph -o control_bin10.bed --binSize 10
sortBed -i control_bin10.bed > control_bin10_sort.bed

### merge columns of bed

### identify p<0.05 peaks using bed 
