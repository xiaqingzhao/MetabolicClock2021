library(data.table) 
library(tidyverse)
library(survival)
library(survminer)
library(scales)
library(splitstackshape) 

DL20<-fread("Lifespan20_Raw.txt")

### Create survival objects 

SurvObj<-function(DL) {
  
  InfoDL<-DL[(Censored>0 | IntDeaths>0), .(AgeH, Chamber, UniqueName, Deaths, Carried, Censored, IntDeaths, N)]
  eDL<-InfoDL[, .(AgeH, Chamber, UniqueName, Censored, IntDeaths)]
  
  Deaths<-expandRows(eDL[,.(AgeH, Chamber, UniqueName, IntDeaths)], "IntDeaths")
  Deaths$Deaths<-1 #Assign value 1 to each expanded row
  Deaths$Censored<-0 #Assign value 0 to each expanded row
  Censored<-expandRows(eDL[,.(AgeH, Chamber, UniqueName, Censored)], "Censored")
  Censored$Censored<-1
  Censored$Deaths<-0
  
  temp<-rbind(Deaths, Censored) #temp$Deaths is equivalent to event in the survival object, as all deaths are 1 and all censored are 0
  Data<-temp[,Genotype:=UniqueName]
  Data$AgeD<-(Data$AgeH)/24
  
  return(Data)
}

#############

ndeath<-rep(NA,length(unique(DL20$UniqueName)))
rmean<-rep(NA,length(unique(DL20$UniqueName)))
se<-rep(NA,length(unique(DL20$UniqueName)))
median<-rep(NA,length(unique(DL20$UniqueName)))

for (i in 1:length(unique(DL20$UniqueName))) {
  
  SubData<-DL20[UniqueName==unique(DL20$UniqueName)[i]]
  Data<-SurvObj(SubData)
  
  fit<-survfit(Surv(AgeD,Deaths)~1, data=Data)
  ndeath[i]<-summary(fit)$table[4]
  rmean[i]<-summary(fit)$table[5]
  se[i]<-summary(fit)$table[6]
  median[i]<-summary(fit)$table[7]
  SummaryStats<-data.table(Genotype=unique(DL20$UniqueName),NDeath=ndeath,Mean=rmean,Se=se,Median=median)
}

SummaryStats<-SummaryStats[order(Genotype)]
write.table(SummaryStats,"MeanLS_20.txt",sep="\t",quote=F,row.names=F)

## Or alternatively, 

Surv20<-SurvObj(DL20)

fit20<-survfit(Surv(AgeD,Deaths)~Genotype,data=Surv20)
print(fit20,rmean="individual")

## Log-rank test comparing survival curves from different genotypes

Diff<-survdiff(Surv(AgeD,Deaths)~Genotype,data=Surv20)
Diff

##############

LS<-SummaryStats
OrderLS<-arrange(LS,Mean)

Color<-seq_gradient_pal("red","blue","Lab")(seq(0,1,length.out=20))
ColorMatch<-data.table(Genotype=OrderLS$Genotype,Color)
ColorMatch<-ColorMatch[order(Genotype)]
OrderedColor<-ColorMatch$Color

fit<-survfit(Surv(AgeD,Deaths)~Genotype, data=Surv20)
p<-ggsurvplot_list(fit,Surv20,xlim=c(0,100),palette=OrderedColor,ggtheme=theme_bw(),legend.title="")+xlab("Time (Days)")+ylab("Survivorship")

pdf("SurvivalCurves_20Genotypes.pdf", width=10, height=4)
p[[1]]+theme(legend.position="none")
dev.off()


pdf("SurvivalCurves_20Genotypes_SampleAges.pdf", width=8, height=4)
p[[1]]+geom_vline(xintercept=c(4,10,24,45,69,80),col="cyan3")+theme(legend.position="none")
dev.off()

######

## Variation across 20 DGRP strains

Data<-SummaryStats
SortData<-Data[order(Mean)]

pdf("LSVariation_20Genotypes.pdf", width=4, height=4)
ggplot(SortData,aes(x=1:20,y=Mean))+geom_point()+geom_errorbar(aes(ymin=Mean-Se,ymax=Mean+Se))+theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  xlab("Genotype")+ylab("Mean Lifespan (Days)")+scale_y_continuous(limits=c(30,85))
dev.off()


######################

## Gompertz-Makeham Parameters vs Lifespan

library(GGally)

Data<-fread("DemographySummary.txt")

Matrix<-as.matrix(Data[,c(3,4,5,6,8,9)])
rownames(Matrix)<-Data$Genotype
colnames(Matrix)<-c("MeanLS","MedianLS","LogAlpha","Beta","LogMiu45","LogMiu60")
pairs(Matrix,pch=20)

DF<-Data[,c(3,5,6)]
a<-ggpairs(DF,axisLabels="none",lower=list(continuous="smooth"),diag=list(continuous="bar"))+theme_bw()

# b<-ggpairs(DF,axisLabels="none",lower=list(continuous="smooth"), columnLabels=c("Mean Lifespan","log (Alpha)","Beta"))+theme_bw()

colnames(DF)<-c("MeanLifespan","'log('*alpha*')'","beta")
b<-ggpairs(DF,axisLabels="none",lower=list(continuous="smooth"), 
           labeller=label_parsed)+theme_bw()
pdf("LSCorrelation.pdf", width=5, height=4)
b
dev.off()

####################

## Which parameter determines mean lifespan? 

lm<-lm(MeanLS~LogMAlpha+MBeta,data=Data)
summary(lm)
anova(lm)

#####################

## Inverse correlation between alpha and beta

cor.test(Data$LogMAlpha,Data$MBeta,method="spearman")

#####################

## Mortality at Day 45 and Day 60 are correlated with LS

cor.test(Data$MeanLS,Data$LogMiu45,method="spearman")
cor.test(Data$MeanLS,Data$LogMiu60,method="spearman")

