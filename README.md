# cuproptosis-and-cervical-cancer
R script code for cuproptosis and cervical cancer.

#1.cuproptosis.cuproptosisExp.R
##if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)            
expFile="symbol.txt"      
geneFile="gene.txt"      
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\08.cuproptosisExp")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

outTab=rbind(ID=colnames(geneExp),geneExp)
write.table(outTab, file="cuproptosisExp.txt", sep="\t", quote=F, col.names=F)

2.cuproptosisLnc.R
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)
corFilter=0.4            
pvalueFilter=0.001       
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\09.cuproptosisLnc")     

rt=read.table("lncRNA.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]
conNum=length(group[group==1])       
treatNum=length(group[group==0])    
sampleType=c(rep(1,conNum), rep(2,treatNum))

rt1=read.table("cuproptosisExp.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
cuproptosis=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
cuproptosis=avereps(cuproptosis)
cuproptosis=cuproptosis[rowMeans(cuproptosis)>0.1,]

group=sapply(strsplit(colnames(cuproptosis),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
cuproptosis=cuproptosis[,group==0]

outTab=data.frame()
for(i in row.names(lncRNA)){
	if(sd(lncRNA[i,])>0.1){
		test=wilcox.test(data[i,] ~ sampleType)
		if(test$p.value<0.05){
			for(j in row.names(cuproptosis)){
				x=as.numeric(lncRNA[i,])
				y=as.numeric(cuproptosis[j,])
				corT=cor.test(x,y)
				cor=corT$estimate
				pvalue=corT$p.value
				if((cor>corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(Cuproptosis=j,lncRNA=i,cor,pvalue,Regulation="postive"))
				}
				if((cor< -corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(Cuproptosis=j,lncRNA=i,cor,pvalue,Regulation="negative"))
				}
			}
		}
	}
}

write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)

cuproptosisLncRNA=unique(as.vector(outTab[,"lncRNA"]))
cuproptosisLncRNAexp=data[cuproptosisLncRNA,]
cuproptosisLncRNAexp=rbind(ID=colnames(cuproptosisLncRNAexp), cuproptosisLncRNAexp)
write.table(cuproptosisLncRNAexp,file="cuproptosisLncExp.txt",sep="\t",quote=F,col.names=F)

3.cuproptosis.ggalluvial.R
install.packages("ggplot2")
install.packages("ggalluvial")

library(dplyr)
library(ggalluvial)
library(ggplot2)
inputFile="corResult.txt"         
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\12.ggalluvial")    

rt=read.table(inputFile, header=T, sep="\t", check.names=F)

mycol=rainbow(length(unique(rt[,"Cuproptosis"])), s=0.8, v=0.8)
#mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)

p1<-ggplot(data=rt, aes(axis1 =lncRNA , axis2 =Cuproptosis, y = 1))+
geom_alluvium(aes(fill =Cuproptosis), width = 0.1, knot.pos = 0.1, reverse = F)+
geom_stratum(fill=NA, color=NA, alpha= 0.5, width = 0.1)+
geom_text(stat = 'stratum', size =1.5, color='black', label.strata = T)+
scale_fill_manual(values = mycol) +
scale_x_discrete(limits = c('lncRNA','Cuproptosis'), expand=c(0, 0))+
xlab("") + ylab("") + theme_bw() +
theme(axis.line = element_blank(), axis.ticks = element_blank(),axis.text.x = element_blank()) +
theme(panel.grid =element_blank()) +
theme(panel.border = element_blank()) +
coord_flip()+ggtitle("")

pdf(file="ggalluvial.pdf", width=8, height=4.75)
print(p1)
dev.off()
4.cuproptosis13.mergeTime.R
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
library(limma)               
lncFile="cuproptosisLncExp.txt"      
cliFile="time.txt"           
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\13.mergeTime")    

rt=read.table(lncFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data=t(data)
data=avereps(data)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.txt",sep="\t",row.names=F,quote=F)

5.cuproptosis14.model.R
install.packages("survival")
install.packages("caret")
install.packages("glmnet")
install.packages("survminer")
install.packages("timeROC")

library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)
coxPfilter=0.05       
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\14.model")      

rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
rt$futime[rt$futime<=0]=1
rt$futime=rt$futime/365
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)

bioForest=function(coxFile=null,forestFile=null,forestCol=null){

rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

pdf(file=forestFile, width=7, height=6)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
LOGindex = 10
hrLow = log(as.numeric(hrLow),LOGindex)
hrHigh = log(as.numeric(hrHigh),LOGindex)
hr = log(as.numeric(hr),LOGindex)
xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), forestCol[1],forestCol[2])
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
a1 = axis(1,labels=F,tick=F)
axis(1,a1,10^a1)
dev.off()
}


n=1000   
for(i in 1:n){

inTrain<-createDataPartition(y=rt[,2], p=0.5, list=F)
train<-rt[inTrain,]
test<-rt[-inTrain,]
trainOut=cbind(id=row.names(train),train)
testOut=cbind(id=row.names(test),test)
#cox
outUniTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(train[,3:ncol(train)])){
#cox
cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]

if(coxP<coxPfilter){
sigGenes=c(sigGenes,i)
outUniTab=rbind(outUniTab,
cbind(id=i,
HR=coxSummary$conf.int[,"exp(coef)"],
HR.95L=coxSummary$conf.int[,"lower .95"],
HR.95H=coxSummary$conf.int[,"upper .95"],
pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
)
}
}
uniSigExp=train[,sigGenes]
uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)
if(ncol(uniSigExp)<6){next}
#lasso
x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
y=data.matrix(Surv(uniSigExp$futime,uniSigExp$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoSigExp=uniSigExp[,c("futime", "fustat", lassoGene)]
lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
if(nrow(geneCoef)<2){next}

multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outMultiTab=data.frame()
outMultiTab=cbind(
coef=multiCoxSum$coefficients[,"coef"],
HR=multiCoxSum$conf.int[,"exp(coef)"],
HR.95L=multiCoxSum$conf.int[,"lower .95"],
HR.95H=multiCoxSum$conf.int[,"upper .95"],
pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
outMultiTab=outMultiTab[,1:2]

riskScore=predict(multiCox,type="risk",newdata=train)      
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))

riskScoreTest=predict(multiCox,type="risk",newdata=test)     #利用train得到模型预测test样品风险
riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest))

diff=survdiff(Surv(futime, fustat) ~risk,data = train)
pValue=1-pchisq(diff$chisq, df=1)
diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
pValueTest=1-pchisq(diffTest$chisq, df=1)
#ROC
predictTime=1    
roc=timeROC(T=train$futime, delta=train$fustat,
marker=riskScore, cause=1,
times=c(predictTime), ROC=TRUE)
rocTest=timeROC(T=test$futime, delta=test$fustat,
marker=riskScoreTest, cause=1,
times=c(predictTime), ROC=TRUE)
if((pValue<0.03) & (roc$AUC[2]>0.65) & (pValueTest<0.05) & (rocTest$AUC[2]>0.63)){

write.table(trainOut,file="data.train.txt",sep="\t",quote=F,row.names=F)
write.table(testOut,file="data.test.txt",sep="\t",quote=F,row.names=F)

write.table(outUniTab,file="uni.trainCox.txt",sep="\t",row.names=F,quote=F)
write.table(uniSigExpOut,file="uni.SigExp.txt",sep="\t",row.names=F,quote=F)
bioForest(coxFile="uni.trainCox.txt",forestFile="uni.foreast.pdf",forestCol=c("red","green"))

write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
pdf("lasso.lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()

write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
write.table(trainRiskOut,file="risk.train.txt",sep="\t",quote=F,row.names=F)
write.table(testRiskOut,file="risk.test.txt",sep="\t",quote=F,row.names=F)

allRiskOut=rbind(trainRiskOut, testRiskOut)
write.table(allRiskOut,file="risk.all.txt",sep="\t",quote=F,row.names=F)
break
}
}

6.cuproptosis.cliStat.R
trainFile="data.train.txt"    #train???????ļ?
testFile="data.test.txt"      #test???????ļ?
cliFile="clinical.txt"        #?ٴ??????ļ?
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\15.cliStat")      #???ù???Ŀ¼


train=read.table(trainFile, header=T, sep="\t", check.names=F, row.names=1)

test=read.table(testFile, header=T, sep="\t", check.names=F, row.names=1)


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))


trainCli=cli[row.names(train),]
trainCli=cbind(trainCli, Type="Train")
testCli=cli[row.names(test),]
testCli=cbind(testCli, Type="Test")
rt=rbind(trainCli, testCli)


cliStatOut=data.frame()
for(i in 1:(ncol(rt)-1)){
	nameStat=colnames(rt)[i]
	tableStat=table(rt[,c(nameStat,"Type")])
	tableStatSum=cbind(Total=rowSums(tableStat), tableStat)
	tableStatRatio=prop.table(tableStatSum,2)
	tableStatRatio=round(tableStatRatio*100,2)
	tableStatPaste=paste(tableStatSum,"(",tableStatRatio,"%)",sep="")
	tableStatOut=matrix(tableStatPaste,ncol=3,dimnames=dimnames(tableStatSum))
	pStat=chisq.test(tableStat[row.names(tableStat)!="unknow",])
	pValueStat=round(pStat$p.value,4)
	pValueCol=c(pValueStat,rep(" ",(nrow(tableStatOut)-1)) )
	tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatOut),tableStatOut,Pvalue=pValueCol)
	cliStatOut=rbind(cliStatOut,tableStatOut)
}


write.table(cliStatOut,file="cliStat.result.xls",sep="\t",quote=F,row.names=F)

7.cuproptosis16.corplot.R
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

install.packages("tidyverse")
install.packages("ggplot2")
install.packages("ggExtra")

library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)

riskFile="risk.all.txt"                
cuprExpFile="cuproptosisExp.txt"       
lncExpFile="cuproptosisLncExp.txt"     
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\16.corplot")     

rt1=read.table(cuprExpFile, header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
cuproptosis=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
cuproptosis=avereps(cuproptosis)
cuproptosis=cuproptosis[rowMeans(cuproptosis)>0.1,]

group=sapply(strsplit(colnames(cuproptosis),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
cuproptosis=t(cuproptosis[,group==0])

rt=read.table(lncExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]

riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
lncRNA=t(lncRNA[colnames(riskRT)[3:(ncol(riskRT)-2)],])

outTab=data.frame()
for(lncrna in colnames(lncRNA)){
	for(gene in colnames(cuproptosis)){
		x=as.numeric(lncRNA[,lncrna])
		y=as.numeric(cuproptosis[,gene])
		corT=cor.test(x, y)
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(cuproptosis=gene, lncrna=lncrna, cor, text, pvalue))
	}
}

outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=7, height=5.6)
ggplot(outTab, aes(lncrna, cuproptosis)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +    
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),   #x??????
	      axis.text.y = element_text(size = 12, face = "bold")) +       #y??????
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #????ͼ??
	scale_x_discrete(position = "bottom")      #????X????????ʾ??λ??
dev.off()

8.cuproptosis.survial.R
install.packages("survival")
install.packages("survminer")

library(survival)
library(survminer)
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\17.survival")     

bioSurvival=function(inputFile=null, outFile=null){
	
	rt=read.table(inputFile, header=T, sep="\t")

	diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           ylab="Overall survival",
		           break.time.by = 1,
		           palette=c("red", "blue"),
		           risk.table=TRUE,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)
	
	pdf(file=outFile, width = 6.5, height =5.5, onefile = FALSE)
	print(surPlot)
	dev.off()
}

bioSurvival(inputFile="risk.train.txt", outFile="surv.train.pdf")
bioSurvival(inputFile="risk.test.txt", outFile="surv.test.pdf")
bioSurvival(inputFile="risk.all.txt", outFile="surv.all.pdf")

9.cuproptosis18.PFS.R
#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

riskFile="risk.all.txt"      
cliFile="Survival_SupplementalTable_S1_20171025_xena_sp"    
setwd("C:\\Users\\lexb\\Desktop\\cuproptosis\\18.PFS")      

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[,c("PFI.time", "PFI")]
cli=na.omit(cli)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365
cli=as.matrix(cli)
row.names(cli)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(cli))

sameSample=intersect(row.names(risk), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], risk[sameSample,"risk",drop=F])

length=length(levels(factor(rt$risk)))
diff=survdiff(Surv(futime, fustat) ~ risk, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit=survfit(Surv(futime, fustat) ~ risk, data = rt)
#print(surv_median(fit))

surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           font.legend=10,
		           xlab="Time(years)",
		           ylab="Progression free survival",
		           break.time.by = 1,
		           palette = c("red", "blue"),
		           risk.table=TRUE,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)

pdf(file="PFS.pdf", width=6.5, height=5.5, onefile=FALSE)
print(surPlot)
dev.off()

10.cuproptosis.riskPlot.R
#install.packages("pheatmap")

library(pheatmap)      
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\19.riskPlot")      


bioRiskPlot=function(inputFile=null, project=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    
	rt=rt[order(rt$riskScore),]      #???ղ??˷??մ??ֶ???Ʒ????????
		
	riskClass=rt[,"risk"]
	lowLength=length(riskClass[riskClass=="low"])
	highLength=length(riskClass[riskClass=="high"])
	lowMax=max(rt$riskScore[riskClass=="low"])
	line=rt[,"riskScore"]
	line[line>10]=10
	pdf(file=paste0(project, ".riskScore.pdf"), width=7, height=4)
	plot(line, type="p", pch=20,
		 xlab="Patients (increasing risk socre)",
		 ylab="Risk score",
		 col=c(rep("blue",lowLength),rep("red",highLength)) )
	abline(h=lowMax,v=lowLength,lty=2)
	legend("topleft", c("High risk","Low Risk"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	dev.off()
		

	color=as.vector(rt$fustat)
	color[color==1]="red"
	color[color==0]="blue"
	pdf(file=paste0(project, ".survStat.pdf"), width=7, height=4)
	plot(rt$futime, pch=19,
		 xlab="Patients (increasing risk socre)",
		 ylab="Survival time (years)",
		 col=color)
	legend("topleft", c("Dead","Alive"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	abline(v=lowLength,lty=2)
	dev.off()
	
	ann_colors=list()
	bioCol=c("blue", "red")
	names(bioCol)=c("low", "high")
	ann_colors[["Risk"]]=bioCol

	rt1=rt[c(3:(ncol(rt)-2))]
	rt1=t(rt1)
	annotation=data.frame(Risk=rt[,ncol(rt)])
	rownames(annotation)=rownames(rt)
	pdf(file=paste0(project, ".heatmap.pdf"), width=7, height=4)
	pheatmap(rt1, 
		     annotation=annotation,
		     annotation_colors = ann_colors, 
		     cluster_cols = FALSE,
		     cluster_rows = FALSE,
		     show_colnames = F,
		     scale="row",
		     color = colorRampPalette(c(rep("blue",3.5), "white", rep("red",3.5)))(50),
		     fontsize_col=3,
		     fontsize=7,
		     fontsize_row=8)
	dev.off()
}

bioRiskPlot(inputFile="risk.train.txt", project="train")

bioRiskPlot(inputFile="risk.test.txt", project="test")

bioRiskPlot(inputFile="risk.all.txt", project="all")

11.cuproptosis20.indep.R
#install.packages('survival')

library(survival)      
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\20.indep")      

bioForest=function(coxFile=null, forestFile=null, forestCol=null){
	#??ȡ?????ļ?
	rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	pdf(file=forestFile, width=6.6, height=4.5)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
		
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=1,col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
	axis(1)
	dev.off()
}

indep=function(riskFile=null, cliFile=null, project=null){
	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #??ȡ?????ļ?
	cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #??ȡ?ٴ??ļ?
	
	
	sameSample=intersect(row.names(cli),row.names(risk))
	risk=risk[sameSample,]
	cli=cli[sameSample,]
	rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
	

	uniCoxFile=paste0(project,".uniCox.txt")
	uniCoxPdf=paste0(project,".uniCox.pdf")
	uniTab=data.frame()
	for(i in colnames(rt[,3:ncol(rt)])){
		 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
		 coxSummary = summary(cox)
		 uniTab=rbind(uniTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
	}
	write.table(uniTab,file=uniCoxFile,sep="\t",row.names=F,quote=F)
	bioForest(coxFile=uniCoxFile, forestFile=uniCoxPdf, forestCol="green")
	
	multiCoxFile=paste0(project,".multiCox.txt")
	multiCoxPdf=paste0(project,".multiCox.pdf")
	uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
	rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
	multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
	multiCoxSum=summary(multiCox)
	multiTab=data.frame()
	multiTab=cbind(
	             HR=multiCoxSum$conf.int[,"exp(coef)"],
	             HR.95L=multiCoxSum$conf.int[,"lower .95"],
	             HR.95H=multiCoxSum$conf.int[,"upper .95"],
	             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	multiTab=cbind(id=row.names(multiTab),multiTab)
	write.table(multiTab, file=multiCoxFile, sep="\t", row.names=F, quote=F)
	bioForest(coxFile=multiCoxFile, forestFile=multiCoxPdf, forestCol="red")
}

indep(riskFile="risk.all.txt", cliFile="clinical.txt", project="all")

12.cuproptosis.ROC.R
#install.packages("survival")
#install.packages("survminer")
install.packages("timeROC")

library(survival)
library(survminer)
library(timeROC)

riskFile="risk.all.txt"     
cliFile="clinical.txt"      
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\21.ROC")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)

ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
	           marker=risk$riskScore,cause=1,
	           weighting='aalen',
	           times=c(1,3,5),ROC=TRUE)
pdf(file="ROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
	   c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	     paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	     paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	   col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()

predictTime=1     
aucText=c()
pdf(file="cliROC.pdf", width=5, height=5)
#???Ʒ??յ÷ֵ?ROC????
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
#???ٴ????ݽ???ѭ?????????ٴ????ݵ?ROC????
for(i in 4:ncol(rt)){
	ROC_rt=timeROC(T=rt$futime,
				   delta=rt$fustat,
				   marker=rt[,i], cause=1,
				   weighting='aalen',
				   times=c(predictTime),ROC=TRUE)
	plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
	aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}

legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()

13.cuproptosis.C-index.R
#install.packages("survival")
#install.packages("rms")
#install.packages("pec")

library(dplyr)
library(survival)
library(rms)
library(pec)

riskFile="risk.all.txt"     
cliFile="clinical.txt"      
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\22.C-index")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)

C-indexֵ
riskScore=cph(Surv(futime,fustat)~riskScore, data=rt, surv=TRUE)
Age=cph(Surv(futime,fustat)~Age, data=rt, surv=TRUE)
Grade=cph(Surv(futime,fustat)~Grade, data=rt, surv=TRUE)
Stage=cph(Surv(futime,fustat)~Stage, data=rt, surv=TRUE)
c_index  <- cindex(list("Risk score"=riskScore, 
                        "Age"=Age,
                        "Grade"=Grade,
                        "Stage"=Stage),
                    formula=Surv(futime,fustat)~ .,
                    data=rt,
                    eval.times=seq(0,10,1),
                    splitMethod="bootcv",
                    B=1000
                    )
pdf(file="C-index.pdf", width=5.5, height=5)
plot(c_index, 
     xlim=c(0,10), ylim=c(0.4,0.8), 
     col=bioCol, xlab="Time (years)",
     legend.x=6, legend.y=0.82, legend.cex=1)
dev.off()

14.cuproptosis23.Nomo.R
#install.packages("survival")
#install.packages("regplot")
#install.packages("rms")

library(survival)
library(regplot)
library(rms)

riskFile="risk.all.txt"      
cliFile="clinical.txt"      
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\23.Nomo")     

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)

res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
              plots = c("density", "boxes"),
              clickable=F,
              title="",
              points=TRUE,
              droplines=TRUE,
              observation=rt[20,],
              rank="sd",
              failtime = c(1,3,5),
              prfail = F)

nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)


pdf(file="calibration.pdf", width=5, height=5)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
	 xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
	   col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

15.cuproptosis24.cliGroupSur.R
#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

riskFile="risk.all.txt"      
cliFile="clinical.txt"       
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\24.cliGroupSur")      

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
			           break.time.by = 2,
			           palette=c("red", "blue"),
			           risk.table=F,
			       	   risk.table.title="",
			           risk.table.col = "strata",
			           risk.table.height=.25)
	
	j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
	pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
			       width = 6,        #ͼƬ?Ŀ???
			       height =5)        #ͼƬ?ĸ߶?
	print(surPlot)
	dev.off()
}

16.cuproptosis.PCA.R
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

install.packages("scatterplot3d")
library(limma)
library(scatterplot3d)
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\25.PCA")       

myPCA=function(input=null,output=null){
	rt=read.table(input, header=T, sep="\t", check.names=F)
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	data=data[rowMeans(data)>0.5,]
	
	type=sapply(strsplit(colnames(data),"\\-"),"[",4)
	type=sapply(strsplit(type,""),"[",1)
	type=gsub("2","1",type)
	data=t(data[,type==0])
	rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
		
	risk=read.table("risk.all.txt", header=T, sep="\t", row.names=1, check.names=F)
	sameSample=intersect(rownames(data),rownames(risk))
	data=data[sameSample,]
	risk=risk[sameSample,]
	group=as.vector(risk[,"risk"])
		
	data.class <- rownames(data)
	data.pca <- prcomp(data, scale. = TRUE)

	color=ifelse(group=="low",4,2)
	pcaPredict=predict(data.pca)
	pdf(file=output, width=7, height=7)
	par(oma=c(1,1,2.5,1))
	s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=35)
	legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c(4,2))
	dev.off()
}

myPCA(input="symbol.txt",output="PCA.allGene.pdf")

myPCA(input="cuproptosisExp.txt",output="PCA.cuproptosisGene.pdf")

myPCA(input="cuproptosisLncExp.txt",output="PCA.cuproptosisLncRNA.pdf")

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
data=risk[,3:(ncol(risk)-2)]
group=as.vector(risk[,"risk"])
		
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)

color=ifelse(group=="low",4,2)
pcaPredict=predict(data.pca)
pdf(file="PCA.riskLnc.pdf", width=6.5, height=6)
par(oma=c(1,1,2.5,1))
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=35)
legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c(4,2))
dev.off()

17.cuproptosis26.riskDiff.R
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)
expFile="symbol.txt"          
riskFile="risk.all.txt"       
logFCfilter=1                 
fdrFilter=0.05                
setwd("C:\\Users\\lexb\\Desktop\\cuproptosis\\26.riskDiff")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(colnames(data), row.names(risk))
data=data[,sameSample]
risk=risk[sameSample,]
riskLow=risk[risk$risk=="low",]
riskHigh=risk[risk$risk=="high",]
dataLow=data[,row.names(riskLow)]
dataHigh=data[,row.names(riskHigh)]
data=cbind(dataLow,dataHigh)
data=data[rowMeans(data)>1,]
conNum=ncol(dataLow)
treatNum=ncol(dataHigh)
Type=c(rep(1,conNum), rep(2,treatNum))

outTab=data.frame()
for(i in row.names(data)){
	rt=data.frame(expression=data[i,], Type=Type)
	wilcoxTest=wilcox.test(expression ~ Type, data=rt)
	pvalue=wilcoxTest$p.value
	conGeneMeans=mean(data[i,1:conNum])
	treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
	logFC=log2(treatGeneMeans)-log2(conGeneMeans)
	conMed=median(data[i,1:conNum])
	treatMed=median(data[i,(conNum+1):ncol(data)])
	diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		outTab=rbind(outTab,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff, file="riskDiff.txt", sep="\t", row.names=F, quote=F)

18.cuproptosis.GO.R
#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("circlize")
#install.packages("RColorBrewer")
#install.packages("ggpubr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("ComplexHeatmap")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library("ggpubr")
library(ComplexHeatmap)

pvalueFilter=0.05       
qvalueFilter=0.05       

#??????ɫ
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\27.GO")      
rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)     

#????????ת??Ϊ????id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

showNum=10
if(nrow(GO)<30){
	showNum=nrow(GO)
}

#??״ͼ
pdf(file="barplot.pdf", width=8, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		
pdf(file="bubble.pdf", width=8, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

data=GO %>% group_by(ONTOLOGY) %>% slice_head(n=10)
pdf(file="barplot.color.pdf", width=8, height=6.5)
ggbarplot(data, x="Description", y="Count", fill = "ONTOLOGY", color = "white",
		  xlab="Term",
          orientation = "horiz",  
          palette = "aaas",    
          legend = "right",    
          sort.val = "asc",    
          sort.by.groups=TRUE)+   
          scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))
dev.off()

ontology.col=c("#00AFBB", "#E7B800", "#90EE90")
data=GO[order(GO$p.adjust),]
datasig=data[data$p.adjust<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

pdf("GO.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()

middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)

main.legend = Legend(
  labels = c("Biological Process","Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()

20.cuproptosis28.KEGG.R
#install.packages("colorspace")
#install.packages("ggplot2")
#install.packages("circlize")
#install.packages("RColorBrewer")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("ComplexHeatmap")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05      
qvalueFilter=1     


colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
setwd("C:\\Users\\z'z\\Desktop\\cuproptosis\\28.KEGG")          
rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)    

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#??״ͼ
pdf(file="barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

pdf(file="bubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()

Pathway.col=c("#90EE90", "#E7B800", "#00AFBB")
showNum=18
data=KEGG[order(KEGG$p.adjust),]
if(nrow(KEGG)>showNum){
	data=data[1:showNum,]
}
data$Pathway="KEGG"
main.col = Pathway.col[as.numeric(as.factor(data$Pathway))]

BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(KEGG=data$ID,start=1,end=max(BgGene))
rownames(df) = df$KEGG
bed2 = data.frame(KEGG=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(KEGG=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(KEGG=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

pdf(file="KEGG.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()

middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',Pathway.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)

main.legend = Legend(
  labels = c("KEGG"),  type="points",pch=15,
  legend_gp = gpar(col=Pathway.col), title_position = "topcenter",
  title = "Pathway", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)

logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()

21.cuproptosis29.immFunction.R
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSEABase")

#install.packages("pheatmap")
#install.packages("reshape2")

library(limma)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(reshape2)

expFile="symbol.txt"         
gmtFile="immune.gmt"         
riskFile="risk.all.txt"      
setwd("C:\\Users\\lexb\\Desktop\\cuproptosis\\29.immFunction")      

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
	
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
	return((x-min(x))/(max(x)-min(x)))}
#??ssGSEA score???н???
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file="immFunScore.txt", sep="\t", quote=F, col.names=F)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
lowSample=row.names(risk[risk$risk=="low",])
highSample=row.names(risk[risk$risk=="high",])
lowData=data[,lowSample]
highData=data[,highSample]
data=cbind(lowData, highData)
conNum=ncol(lowData)        #?ͷ???????Ʒ??Ŀ
treatNum=ncol(highData)     #?߷???????Ʒ??Ŀ
sampleType=c(rep(1,conNum), rep(2,treatNum))

sigVec=c()
for(i in row.names(data)){
	test=wilcox.test(data[i,] ~ sampleType)
	pvalue=test$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	sigVec=c(sigVec, paste0(i, Sig))
}
row.names(data)=sigVec

Type=c(rep("Low risk",conNum), rep("High risk",treatNum))
Type=factor(Type, levels=c("Low risk", "High risk"))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf("heatmap.pdf", width=8, height=4.6)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()

22.cuproptosis.maftools.R
#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")

library(maftools)      
setwd("C:\\Users\\lexb\\Desktop\\cuproptosis\\31.maftools")     

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F)
outTab=risk[,c(1, ncol(risk))]
colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)

geneNum=15     
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

#ע?͵???ɫ
ann_colors=list()
col=c("blue", "red")
names(col)=c("low", "high")
ann_colors[["Risk"]]=col

pdf(file="low.pdf", width=6, height=6)
maf=read.maf(maf="low.maf", clinicalData="ann.txt")    
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()


pdf(file="high.pdf", width=6, height=6)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")    
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

23.cuproptosis.riskTMB.R
#install.packages("ggpubr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)
library(ggpubr)
setwd("C:\\Users\\lexb\\Desktop\\cuproptosis\\32.riskTMB")     

tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)
	
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
	
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)
data$TMB=log2(data$TMB+1)
	
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
boxplot=ggviolin(data, x="risk", y="TMB", fill="risk",
			      xlab="",
			      ylab="Tumor tmbation burden (log2)",
			      legend.title="",
			      palette = c("#0066FF","#FF0000"),
			      add = "boxplot", add.params = list(fill="white"))+ 
	stat_compare_means(comparisons = my_comparisons)
	
pdf(file="riskTMB.pdf", width=5, height=4.5)
print(boxplot)
dev.off()

24.cuproptosis.tmbSur.R
#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

tmbFile="TMB.txt"            #????ͻ?为???ļ?
riskFile="risk.all.txt"      #?????ļ?
setwd("C:\\Users\\lexb\\Desktop\\cuproptosis\\33.tmbSur")      #???ù???Ŀ¼


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #??ȡ?????ļ?
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)      #??ȡ????ͻ?为???ļ?


sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, tmb)

res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
scoreType=ifelse(data$risk=="low", "low risk", "high risk")
mergeType=paste0(tmbType, "+", scoreType)

bioSurvival=function(surData=null, outFile=null){
   
	diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
	length=length(levels(factor(surData[,"group"])))
	pValue=1-pchisq(diff$chisq, df=length-1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
	#print(surv_median(fit))
	

	bioCol=c("#FF0000","#0066FF","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
	bioCol=bioCol[1:length]
	surPlot=ggsurvplot(fit, 
			           data=surData,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           legend.title="",
			           legend.labs=levels(factor(surData[,"group"])),
			           font.legend=10,
			           legend = c(0.8, 0.8),
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette = bioCol,
			           surv.median.line = "hv",
			           risk.table=F,
			           cumevents=F,
			           risk.table.height=.25)
		pdf(file=outFile, onefile = FALSE, width=5.5, height=4.8)
	print(surPlot)
	dev.off()
}

data$group=tmbType
bioSurvival(surData=data, outFile="TMB.survival.pdf")

data$group=mergeType
bioSurvival(surData=data, outFile="TMB-risk.survival.pdf")

25.cuproptosis.TIDE.R
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

library(limma)
library(ggpubr)
tideFile="TIDE.txt"         
riskFile="risk.all.txt"     
setwd("C:\\Users\\lexb\\Desktop\\34.TIDE")    

tide=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0,,drop=F]
tide=as.matrix(tide)
row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
tide=avereps(tide)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	
sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(tide, risk)
	
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

gg1=ggviolin(data, x="risk", y="TIDE", fill = "risk", 
	         xlab="", ylab="TIDE",
	         palette=c("#0066FF","#FF0000"),
	         legend.title="Risk",
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

pdf(file="TIDE.pdf", width=6, height=5)
print(gg1)
dev.off()

26.cuproptosis.pRRophetic.R
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("car", "ridge", "preprocessCore", "genefilter", "sva"))

#install.packages("ggpubr")

library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)

pFilter=0.001               #pvalue?Ĺ???????
expFile="symbol.txt"        #?????????ļ?
riskFile="risk.all.txt"     #?????ļ?
setwd("C:\\Users\\lexb\\Desktop\\cuproptosis\\35.pRRophetic")     #???ù???Ŀ¼

data(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
allDrugs=unique(drugData2016$Drug.name)

rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)

riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
riskRT$riskScore[riskRT$riskScore>quantile(riskRT$riskScore,0.99)]=quantile(riskRT$riskScore,0.99)


for(drug in allDrugs){
	#Ԥ??ҩ????????
	possibleError=tryCatch(
    	{senstivity=pRRopheticPredict(data, drug, selection=1, dataset = "cgp2016")},
    	error=function(e) e)
    if(inherits(possibleError, "error")){next}
	senstivity=senstivity[senstivity!="NaN"]
	senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
	
	sameSample=intersect(row.names(riskRT), names(senstivity))
	risk=riskRT[sameSample, c("riskScore","risk"),drop=F]
	senstivity=senstivity[sameSample]
	rt=cbind(risk, senstivity)
	
	rt$risk=factor(rt$risk, levels=c("low", "high"))
	type=levels(factor(rt[,"risk"]))
	comp=combn(type, 2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	

	test=wilcox.test(senstivity~risk, data=rt)
	diffPvalue=test$p.value

	x=as.numeric(rt[,"riskScore"])
	y=as.numeric(rt[,"senstivity"])
	corT=cor.test(x, y, method="spearman")
	corPvalue=corT$p.value
	
	if((diffPvalue<pFilter) & (corPvalue<pFilter)){
	
		boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",
					      xlab="Risk",
					      ylab=paste0(drug, " senstivity (IC50)"),
					      legend.title="Risk",
					      palette=c("#0066FF","#FF0000")
					     )+ 
			stat_compare_means(comparisons=my_comparisons)
		pdf(file=paste0("durgSenstivity.", drug, ".pdf"), width=5, height=4.5)
		print(boxplot)
		dev.off()
	
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
				xlab("Risk score") + ylab(paste0(drug, " senstivity (IC50)"))+
				geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				stat_cor(method = 'spearman', aes(x =x, y =y))
	
		pdf(file=paste0("Cor.", drug, ".pdf"), width=5, height=4.6)
		print(p1)
		dev.off()
	}
}

27.cuproOmics.vioplot.R
#install.packages("reshape2")
#install.packages("ggpubr")

library(reshape2)
library(ggpubr)

riskFile="risk.all.txt"      
TMEfile="TMEscores.txt"      
setwd("C:\\biowolf\\cuproOmics\\45.estimateVioplot")     

Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
Risk$risk=factor(Risk$risk, levels=c("low","high"))

score=read.table(TMEfile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:3]
rownames(score)=gsub("(.*?)\\_(.*?)", "\\2", rownames(score))
score=score[row.names(Risk),,drop=F]

#???ݺϲ?
rt=cbind(Risk[,"risk",drop=F], score)

data=melt(rt, id.vars=c("risk"))
colnames(data)=c("Risk", "scoreType", "Score")

p=ggviolin(data, x="scoreType", y="Score", fill = "Risk",
	     xlab="",
	     ylab="TME score",
	     legend.title="Risk",
	     add = "boxplot", add.params = list(color="white"),
	     palette = c("#0088FF", "#FF5555"), width=1)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Risk),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

pdf(file="vioplot.pdf", width=6, height=5)
print(p1)
dev.off()


















































