# R script for GLM analysis of Light/Dark RNA pol II ChIPseq
library("edgeR")
stamp <- format(Sys.time(), "%m-%d-%y");

x <- read.csv("count.csv",row.names="Gene")
group <- factor(c(2,2,2,3,3,3))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)

design <- model.matrix(~group)

y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

fit <- glmFit(y,design)
lrt.light <- glmLRT(fit,coef=2)
topTags(lrt.light)
lrt.dark <- glmLRT(fit,coef=3)
topTags(lrt.dark)

# Dark vs Light. Negative logFC indicates up-regulation in light
lrt.lvsd <- glmLRT(fit,contrast=c(0,-1,1))
tt <- topTags(lrt.lvsd, n=dim(y)[1], sort.by="PValue")
write.table(tt, file=paste("dark-v-light", stamp, "txt", sep=".") , sep="\t")
