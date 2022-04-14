#install.packages("survival")
#install.packages("survminer")

setwd("H:\\000Bio\\TCGA_liver cancer\\12.survival\\STIL")             
library(survival)
library("survminer")
rt=read.table("STIL_TCGA.txt",header=T,sep="\t")                     
#rt=read.table("STIL-GSE.txt",header=T,sep="\t")                       
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)


#survival plot
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(rt[,"risk"])))]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.labs=levels(factor(rt[,"risk"])),
                   legend.title="risk",
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="STIL-GSE.os.pdf", onefile = FALSE, width=8, height=6.5)
print(surPlot)
dev.off()

summary(fit)    
