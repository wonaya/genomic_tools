input <- read.delim("VE_sequencability_normalized_signal.bedgraph", sep="\t") 
bins=51

colnames(input)<-c("chr","start","end","input")
my.half <- floor(bins/2)
my.mid <- ceiling(bins/2)
my.min.max <- matrix(nrow=nrow(input), ncol=1)
colnames(my.min.max) <- c("input_max")
for(n in my.half:nrow(input))
    {
    my.min.max[n,1] <- which.max(input$input[(n-my.half): (n+my.half)])==my.mid
    }
input_min.max <- cbind(input, my.min.max)
input_max <- subset(input_min.max, input_max==TRUE, select=c("chr","start","end","input"))
write.table(input_max, file="VE_max_2.bed", sep="\t", row.names = FALSE, col.names= FALSE)
