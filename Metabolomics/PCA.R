library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)

load("CenterScaledData.RData")
SampleInfo<-fread("Metabolomics1_Randomized.txt")
SampleInfo<-SampleInfo[SampleID!="Sample181",]
SampleInfo[Age=="Day04", AgeNum:=4]
SampleInfo[Age=="Day10", AgeNum:=10]
SampleInfo[Age=="Day24", AgeNum:=24]
SampleInfo[Age=="Day45", AgeNum:=45]
SampleInfo[Age=="Day69", AgeNum:=69]
SampleInfo[Age=="Day80", AgeNum:=80]

## Create PCA plot for publication, using CenterScaleBySample

PCA<-prcomp(CenterScaleBySample,scale.=TRUE)
DF<-data.frame(Genotype=SampleInfo$Genotype,Age=SampleInfo$Age,
               Vial=SampleInfo$VialNum,ID=SampleInfo$SampleID,PCA$x[,1:3])

percentage<-round(PCA$sdev^2/sum(PCA$sdev^2)*100, 2)
percentage<-paste(colnames(PCA$x), "(", as.character(percentage), "%)",sep="" )

pdf("PCA_Age_AllData_CenterScaleBySample.pdf",width=8,height=3)
ggplot(data=DF, aes(x=PC1, y=PC2,colour=Age))+geom_point()+
  xlab(percentage[1])+ylab(percentage[2])+
  theme_bw()+stat_ellipse(level=0.9)+
  scale_color_manual(values=c("red", "orange", "khaki","green","blue","purple"))
dev.off()

####

## PC1 vs Age 

AgeNum<-as.numeric(str_remove(DF$Age,"Day"))
DF<-cbind(DF,AgeNum)

pdf("PC1vsAge.pdf",width=3,height=3)
ggplot(data=DF,aes(x=AgeNum,y=PC1))+geom_point(size=0.6)+theme_bw()+
  xlab("Sample Age")+ylab("PC1")
dev.off()

########

## Lifespan Groups

LS<-fread("LifeSpan.txt")
LS[Mean<50,Group:="Short"]
LS[Mean>50 & Mean<65,Group:="Medium"]
LS[Mean>65,Group:="Long"]
LSGroup<-LS[,c(1,7)]

SampleInfo<-setDT(LSGroup)[SampleInfo,on="Genotype"]

library(plotrix)

DF<-data.frame(SampleInfo,PCA$x[,1:6])

PC1_MeanSE<-DF %>% group_by(Group,Age) %>% summarise(Mean=mean(PC1),SE=std.error(PC1))
PC1_MeanSE<-data.frame(PC1_MeanSE)

pdf("PC1_LongMediumShortGenotype_SE.pdf",width=5,height=4)
ggplot(PC1_MeanSE,aes(x=Age,y=Mean,group=Group,color=Group))+geom_line()+geom_point(position=position_dodge(0.15))+
  geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),width=0.2,position=position_dodge(0.15))+
  theme_bw()+scale_color_manual(values=c('blue','#999999','red'))+
  labs(y="PC1")+guides(col=guide_legend("Lifespan"))
dev.off()

#########################

## Linear model on PC1

fit<-lm(PC1~AgeNum*Group,data=DF)
summary(fit)
anova(fit)


