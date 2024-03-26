args = commandArgs(trailingOnly=TRUE)

input_file<-args[1]
output_file_depth<-args[2]
output_file_flagstats<-args[3]

input_table<-read.table(input_file, sep="\t", header=TRUE)
pdf(output_file_depth, width=8, height=6)
boxplot(input_table$raw_mean_depth, input_table$clean_mean_depth, names =c("raw","cleaned"))
title("Mean sequencing depth (genome-wide)")
stripchart(list(input_table$raw_mean_depth, input_table$clean_mean_depth), vertical=TRUE, add=TRUE, method="jitter", pch=19, col=4, cex=0.8)
dev.off()
pdf(output_file_flagstats, width=8, height=6)
par(mar=c(8.1, 4.1, 4.1, 2.1))
boxplot(input_table$total, input_table$mapped, input_table$duplicates, input_table$paired, input_table$properly_paired, input_table$mate_diff_chr, input_table$mate_diff_chr_mq5, names=c("All reads", "Mapped reads", "Duplicates", "Paired", "Properly paired", "Mate on\ndifferent Chr.", "Mate on different\nChr., MQ>=5" ), las=2)
title("Flagstats")
stripchart(list(input_table$total, input_table$mapped, input_table$duplicates, input_table$paired, input_table$properly_paired, input_table$mate_diff_chr, input_table$mate_diff_chr_mq5), vertical=TRUE, add=TRUE, method="jitter", pch=19, col=4, cex=0.6)
dev.off()
