#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")
#install.packages("reshape2")



library(limma)
library(ggpubr)
library(reshape2)

riskFile="risk.all.txt"         
scoreFile="ssgseaOut.txt"       
setwd("C:\\23.ssGSEAdiff")      


data=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0,drop=F]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data=avereps(t(data))


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,]
risk=risk[sameSample,]
rt=cbind(data,risk[,c("riskScore","risk")])
rt=rt[,-(ncol(rt)-1)]


immCell=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages",
          "Mast_cells","Neutrophils","NK_cells","pDCs","T_helper_cells",
          "Tfh","Th1_cells","Th2_cells","TIL","Treg")
rt1=rt[,c(immCell,"risk")]
data=melt(rt1,id.vars=c("risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))
p=ggboxplot(data, x="Type", y="Score", color = "Risk",
     xlab="",ylab="Score",add = "none",palette = c("blue","red") )
p=p+rotate_x_text(50)

pdf(file="immCell.boxplot.pdf",width=7,height=6)
p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()


immFunction=c("APC_co_inhibition","APC_co_stimulation","CCR",
          "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
          "MHC_class_I","Parainflammation","T_cell_co-inhibition",
          "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
rt1=rt[,c(immFunction,"risk")]
data=melt(rt1,id.vars=c("risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))
p=ggboxplot(data, x="Type", y="Score", color = "Risk",
     xlab="",ylab="Score",add = "none",palette = c("blue","red") )
p=p+rotate_x_text(50)

pdf(file="immFunction.boxplot.pdf",width=7,height=6)
p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()
