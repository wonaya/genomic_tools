library(bsseq)
x <- read.bismark(files = c("CpG_context_B73_ACAGTG_merged_chr10.bedGraph", "CpG_context_Mo17_GCCAAT_merged_chr10.bedGraph"), sampleNames= c("B73", "Mo17"), rmZeroCov=FALSE, verbose=TRUE)
x
seqnames(x)
sample <- BSmooth(x, parallelBy="sample", mc.cores=12, verbose=TRUE)
data(sample)
sample
