#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")



library(limma)
library(pheatmap)

expFile="ARGLncExp.txt"           
uniCoxFile="uni.trainCox.txt"     
setwd("C:\\10.uniSigHeatmap")      


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)


uniCox=read.table(uniCoxFile, header=T, sep="\t", check.names=F, row.names=1)
data=data[row.names(uniCox),]


group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       
treatNum=length(group[group==0])     
sampleType=c(rep(1,conNum), rep(2,treatNum))


sigVec=c()
outTab=data.frame()
for(i in rownames(data)){
	if(sd(data[i,])<0.001){next}
	wilcoxTest=wilcox.test(data[i,] ~ sampleType)
	pvalue=wilcoxTest$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	sigVec=c(sigVec, paste0(i, Sig))
}


exp=log2(data+0.1)
row.names(exp)=sigVec
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=9, height=6)
pheatmap(exp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames = F,
         show_rownames = T,
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()
