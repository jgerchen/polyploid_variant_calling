args = commandArgs(trailingOnly=TRUE)
args_input_table<-read.table(args[1], na.strings = c("NA", "."))

args_MQ_threshold=as.numeric(args[2])
args_MQ_plot=args[3]

args_QD_threshold<-as.numeric(args[4])
args_QD_plot=args[5]

args_FS_threshold<-as.numeric(args[6])
args_FS_plot<-args[7]

args_MQRankSum_threshold<-as.numeric(args[8])
args_MQRankSum_plot<-args[9]

args_ReadPosRankSum_threshold<-as.numeric(args[10])
args_ReadPosRankSum_plot<-args[11]

args_SOR_threshold<-as.numeric(args[12])
args_SOR_plot<-args[13]
#plot MQ
pdf(args_MQ_plot, 10, 7)
hist(args_input_table$V4, xlab="QD", main="RMS Mapping Quality", nclass=100)
abline(v=args_MQ_threshold, col="red", lwd=2, lty=3)
dev.off()
#plot QD
pdf(args_QD_plot, 10, 7)
hist(args_input_table$V5, xlab="QD", main="Variant Confidence/Quality by Depth", nclass=100)
abline(v=args_QD_threshold, col="red", lwd=2, lty=3)
dev.off()
#plot FS
pdf(args_FS_plot, 10, 7)
hist(args_input_table$V6, xlab="FS", main="Phred-scaled p-value using Fisher's exact test to detect strand bias", nclass=100)
abline(v=args_FS_threshold, col="red", lwd=2, lty=3)
dev.off()
#plot MQRankSum
pdf(args_MQRankSum_plot, 10, 7)
hist(args_input_table$V7, xlab="MQRankSum", main="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities" , nclass=100)
abline(v=args_MQRankSum_threshold, col="red", lwd=2, lty=3)
dev.off()
#plot ReadPosRankSum
pdf(args_ReadPosRankSum_plot, 10, 7)
hist(args_input_table$V8, xlab="ReadPosRankSum", main="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias" , nclass=100)
abline(v=args_ReadPosRankSum_threshold, col="red", lwd=2, lty=3)
dev.off()
#plot SOR
pdf(args_SOR_plot, 10, 7)
hist(args_input_table$V9, xlab="SOR", main="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias" , nclass=100)
abline(v=args_SOR_threshold, col="red", lwd=2, lty=3)
dev.off()
