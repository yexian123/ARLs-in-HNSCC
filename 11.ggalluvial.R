#install.packages("ggplot2")
#install.packages("ggalluvial")



library(ggalluvial)
library(ggplot2)
library(dplyr)

corFile="net.network.txt"         
multiCoxFile="multiCox.txt"     
setwd("C:\\11.ggalluvial")     


rt=read.table(corFile, header=T, sep="\t", check.names=F)
uniCox=read.table(multiCoxFile, header=T, sep="\t", check.names=F, row.names=1)


rt=rt[which(rt[,2] %in% row.names(uniCox)),]
rt=rt[,c("ARG","lncRNA","Regulation")]
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")


pdf(file="ggalluvial.pdf", width=7, height=6)
mycol=rep(c("#223D6C","#D8D155","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767","#D20A13","#431A3D","#91612D","#FFD121","#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 
  	 geom_flow(width = 2/10,aes.flow = "forward") + 
	 geom_stratum(alpha = .9, width = 3/10) +
	 scale_fill_manual(values = mycol) +
	 
	 geom_text(stat = "stratum", size = 3, color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #È¥µô×ø±êÖá
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = "none")                            
dev.off()
