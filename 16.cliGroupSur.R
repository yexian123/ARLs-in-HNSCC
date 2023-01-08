#install.packages("survival")
#install.packages("survminer")



library(survival)
library(survminer)

riskFile="risk.all.txt"      
cliFile="N.txt"       
setwd("C:\\16.cliGroupSur")       


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]


sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]


for(j in names(tab)){
	rt1=rt[(rt[,"clinical"]==j),]
	tab1=table(rt1[,"Risk"])
	tab1=tab1[tab1!=0]
	labels=names(tab1)
	if(length(labels)!=2){next}
	if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
		titleName=paste0("age",j)
	}
	
	
	diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	
	
	fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
	surPlot=ggsurvplot(fit, 
			           data=rt1,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           title=paste0("Patients with ",j),
			           legend.title="Risk",
			           legend.labs=labels,
			           font.legend=12,
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette=c("red", "blue"),
			           risk.table=F,
			       	   risk.table.title="",
			           risk.table.col = "strata",
			           risk.table.height=.25)
	
	
	j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
	pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
			       width = 6,        
			       height =5)        
	print(surPlot)
	dev.off()
}
