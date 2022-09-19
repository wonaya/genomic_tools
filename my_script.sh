#!/bin/bash
#SBATCH -J CHH_5       # Job Name
#SBATCH -o CHH_5.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 1             # Total number of mpi tasks requested
#SBATCH -p largemem       # Queue (partition) name -- normal, development, etc.
#SBATCH -t 12:00:00      # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A Springer_Vaughn       # Project_name

#### INPUT SCRIPT HERE ###
#/scratch/02114/wonaya/PeakRanger/peakranger ccat --control Samples1_bwa_unique.bam --format bam --data Samples3_bwa_unique.bam --output Samples3_ccat_result --verbose
#module load bowtie
#module load bismark
module unload R
module load R_mkl
source ~/.bashrc
Rscript r_test.R CHH_context_B73_ACAGTG_merged_chr5.methylKit CHH_context_Mo17_GCCAAT_merged_chr5.methylKit CHH_context_chr5.out
#date
#python extract_regions.py
#date
#python ontology.py Zm_B73_5b_FGS_cds_2012.txt
#bismark -n 1 /scratch/02114/wonaya/TAIR10 -1 LID57018_2_ATGAGC_L007_R1_001.fastq -2 LID57018_2_ATGAGC_L007_R2_001.fastq -p 2 --bowtie2
