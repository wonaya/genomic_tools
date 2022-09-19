#!/bin/bash
 
for C in *_S*-bwa-pe.bam
do
BASE=$(basename $C .bam)
echo "Processing $BASE..."
samtools sort -n -m 1000000000 $C ${BASE}.namesorted
samtools view ${BASE}.namesorted.bam | htseq-count --mode=union --stranded=no -a 10 -t "protein_coding_gene" -i "ID" - Arabidopsis_thaliana.TAIR10.19.protein_coding_gene.gff3 > ${BASE}.namesorted.counts
echo "Done."
done