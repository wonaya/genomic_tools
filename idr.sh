## for individual
samtools view -b -F 1548 -q 30 chipSampleRep1.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > chipSampleRep1.tagAlign.gz
## pool control into one
zcat controlSampleRep1.tagAlign.gz controlSampleRep2.tagAlign.gz controlSampleRep3.tagAlign.gz | gzip -c > controlSampleRep0.tagAlign.gz
## macs2 on individual
macs2 callpeak -t  2_S2_L001_val.sam.sorted.tagAlign.gz -c 1_S1_L001_val.sam.sorted.tagAlign.gz -f BED -n 2_S2_vs_1_S1 -g hs -p 1e-3 --to-large --nomodel --shiftsize 75
macs2 callpeak -t  7_S7_L001_val.sam.sorted.tagAlign.gz -c 1_S1_L001_val.sam.sorted.tagAlign.gz -f BED -n 7_S7_vs_1_S1 -g hs -p 1e-3 --to-large --nomodel --shiftsize 75
macs2 callpeak -t  12_S12_L001_val.sam.sorted.tagAlign.gz -c 1_S1_L001_val.sam.sorted.tagAlign.gz -f BED -n 12_S12_vs_1_S1 -g hs -p 1e-3 --to-large --nomodel --shiftsize 75
## sort peaks
sort -k 8nr,8nr 2_S2_vs_1_S1_peaks.narrowPeak | head -n 100000 | gzip -c > 2_S2_vs_1_S1.regionPeak.gz
sort -k 8nr,8nr 7_S7_vs_1_S1_peaks.narrowPeak | head -n 100000 | gzip -c > 7_S7_vs_1_S1.regionPeak.gz
sort -k 8nr,8nr 12_S12_vs_1_S1_peaks.narrowPeak | head -n 100000 | gzip -c > 12_S12_vs_1_S1.regionPeak.gz

## pool samples
zcat 2_S2_L001_val.sam.sorted.tagAlign.gz 7_S7_L001_val.sam.sorted.tagAlign.gz 12_S12_L001_val.sam.sorted.tagAlign.gz | gzip -c > 2_7_12.tagAlign.gz
## call on pooled sample
macs2 callpeak -t 2_7_12.tagAlign.gz -c 1_S1_L001_val.sam.sorted.tagAlign.gz -f BED -n 2_7_12_vs_1_S1 -g hs -p 1e-3 --to-large --nomodel --shiftsize 75

## between replicates
Rscript batch-consistency-analysis.r /peaks/reps/chipSampleRep1_VS_controlSampleRep0.regionPeak /peaks/reps/chipSampleRep2_VS_controlSampleRep0.regionPeak -1 /consistency/reps/chipSampleRep1_VS_chipSampleRep2 0 F p.value