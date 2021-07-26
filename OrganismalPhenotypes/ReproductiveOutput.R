library(data.table)
library(tidyverse)

Fec<-fread("ReproductiveOutput_Raw.txt")
Fec[,PCF:=Fecundity/Interval/N]
Fec[,ReproductiveOutput:=sqrt(PCF)]

Fec<-Fec[Time!="Day04"]

shapiro.test(Fec$ReproductiveOutput)

Fec$Vial<-as.character(rep(1:99,2))

Fec08<-filter(Fec,Time=="Day08")
Fec12<-filter(Fec,Time=="Day12")

#######

anova(lm(ReproductiveOutput~Block+Genotype+Time+Vial,data=Fec))
## So vials contributes little variance
## It's OK if we take the average of the five replicate vials

#######

by_genotype<-group_by(Fec,Genotype,Time,Block) 
SumFec<-summarise(by_genotype,VialN=n(),MeanRO=mean(ReproductiveOutput,na.rm=TRUE),SdRO=sd(ReproductiveOutput,na.rm=TRUE))
SumFec<-data.table(SumFec)

## Substantial genetic variation 
anova(lm(MeanRO~Block+Genotype+Time,data=SumFec))

Short<-select(SumFec,Genotype,Block,Time,MeanRO)
Short<-spread(Short,Time,MeanRO)
colnames(Short)[3:4]<-c("RepOutput_Day08","RepOutput_Day12")
write.table(Short,"ReproductiveOutputInfo_Genotype20.txt",quote=F,sep="\t",row.names=F)

## Correlation between Day 8 and Day 12 reproductive output
cor.test(Short$RepOutput_Day08, Short$RepOutput_Day12)

############

## No block effects

fit<-lm(RepOutput_Day08~Block,data=Short)
summary(fit)

fit<-lm(RepOutput_Day12~Block,data=Short)
summary(fit)

## So experiment block doesn't contribute to much variation 

###########

## Plot to show variation 

Data<-Short
Data<-Data[,-2]
setDT(Data)[,TotalRO:=RepOutput_Day08+RepOutput_Day12]
Data_Ordered<-Data[order(TotalRO)]
setDT(Data_Ordered)[,GenotypeID:=sprintf("%03d",seq(1:20))]

colnames(Data_Ordered)[2:3]<-c("Day08","Day12")
Long<-gather(Data_Ordered, Age, RO, Day08:Day12)

pdf("ReproductiveOutput_Genotype20.pdf", width=8, height=4)
y_title<-expression(paste(italic("per capita")," Reproductive Output"))
ggplot(Long,aes(fill=factor(Age,levels=c("Day12","Day08")),y=RO,x=GenotypeID))+
  geom_bar(position="stack",stat="identity")+theme_bw()+
  theme(legend.position=c(0.1,0.8),axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Genotype")+
  ylab(y_title)+labs(fill="Age")
dev.off()

pdf("ReproductiveOutputCorrelation_Genotype20.pdf", width=3, height=3)
ggplot(Data,aes(x=RepOutput_Day08,y=RepOutput_Day12))+geom_point()+theme_bw()+
  xlab("Reproductive Output at Day 8")+ylab("Reproductive Output at Day 12")
dev.off()

############

## Day 08 reproductive output correlates with Day 12 reproductive output

cor.test(Data$RepOutput_Day08,Data$RepOutput_Day12)

