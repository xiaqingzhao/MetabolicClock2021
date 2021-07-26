library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(gridExtra)
library(plotly)

Data<-fread("TargetedMetabolomicsAug2019.csv")
Data<-rename(Data,Metabolite=Molecule,Formula='Molecule Formula',SampleID='Replicate Name',
             NormMethod='Normalization Method',RT='Average Measured Retention Time',
             TotalArea='Total Area Fragment',NormArea='Normalized Area')
Data[NormMethod=="Ratio to Global Standards", NormMethod:="GlobalStandards"]
Data[NormMethod=="Ratio to 13C", NormMethod:="RatioTo13C"]
Data<-Data[SampleID!="Sample181",]

SampleInfo<-fread("Metabolomics1_Randomized.txt")
SampleInfo<-SampleInfo[SampleID!="Sample181",]

SampleID<-unique(Data$SampleID)
MetID<-unique(Data$Metabolite)
Matrix<-matrix(Data$TotalArea,nrow=196) 
rownames(Matrix)<-SampleID
colnames(Matrix)<-MetID

Raw<-Matrix[1:181,]

LogRaw<-log10(Raw)

SampleAve<-apply(LogRaw,1,mean)
mean(SampleAve)
sd(SampleAve)
plot(density(SampleAve))

CenterAfterLog<-LogRaw-SampleAve

## Center and scale by sample 

SampleSD<-apply(LogRaw,1,sd)
CenterScaleBySample<-(LogRaw-SampleAve)/SampleSD

save(CenterScaleBySample,file="CenterScaledData.RData")


