#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")



library(limma)
library(ggpubr)
scoreFile="TMEscores.txt"     
riskFile="risk.all.txt"       
setwd("C:\\25.TMEdiff")     


rt=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(rt)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,"risk",drop=F]
rt=cbind(data, risk)


rt$risk=factor(rt$risk, levels=c("low", "high"))
group=levels(factor(rt$risk))
rt$risk=factor(rt$risk, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


for(i in colnames(rt)[1:3]){
	boxplot=ggboxplot(rt, x="risk", y=i, fill="risk",
			          xlab="",
			          ylab=i,
			          legend.title="Risk",
			          palette=c("#0066FF","#FF0000")
			          )+ 
		    stat_compare_means(comparisons=my_comparisons)

	
	pdf(file=paste0(i, ".pdf"), width=5, height=4.5)
	print(boxplot)
	dev.off()
}

