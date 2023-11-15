
args = commandArgs(trailingOnly=TRUE)
sample_name<-args[2]
input_file<-args[1]
output_file<-args[3]

input_flagstats<-read.table(file = input_file, header=FALSE, sep="\t", as.is = TRUE)
row_to_include<-c(1,2,4,5,6,7,9,11,12,13,14,16,17,19,20)

pdf(output_file)
bp<-barplot(as.numeric(input_flagstats$V1[row_to_include]), horiz = TRUE, width=0.1, border=NA)
text(x=max(as.numeric(input_flagstats$V1[row_to_include]))/2, y=bp, paste(input_flagstats$V3[row_to_include], input_flagstats$V1[row_to_include]))
title(sample_name)
dev.off()
