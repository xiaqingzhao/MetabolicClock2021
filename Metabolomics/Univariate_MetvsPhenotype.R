library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(corrplot)
library(circlize)
library(ComplexHeatmap)
library(lineup)

load("CenterScaledData.RData")
Matrix<-CenterScaleBySample

SampleInfo<-fread("Metabolomics1_Randomized.txt")
SampleInfo<-select(SampleInfo,SampleID,Genotype,Age,MetBatch,ExpBlock,VialCopyID)
SampleInfo<-SampleInfo[-181,]

Pheno<-fread("Pheno20.txt")
Pheno<-select(Pheno,Genotype,Block,MeanLS,LogMAlpha,MBeta)
colnames(Pheno)[4:5]<-c("logAlpha","Beta")

PhenoMatrix<-select(Pheno,-(Genotype:Block))
PhenoMatrix<-as.matrix(PhenoMatrix)
rownames(PhenoMatrix)<-Pheno$Genotype

LS<-fread("Pheno20.txt")
LS<-select(LS,Genotype,Block,MeanLS,LogMAlpha,MBeta)
colnames(LS)[4:5]<-c("logAlpha","Beta")
SampleInfo.LS<-merge(SampleInfo,LS,by='Genotype')
SampleInfo.LS<-SampleInfo.LS[match(rownames(Matrix),SampleID),]

##################################################

# Simple linear model of metabolite levels at individule ages vs. LS

LM_Met_LS<-function(Matrix,Phenotype,A) {
  Subset<-Matrix[SampleInfo$Age==A,]
  SampleInfo.LS.Subset<-SampleInfo.LS[SampleInfo$Age==A,]
  Pheno<-SampleInfo.LS.Subset[,..Phenotype]
  Pheno<-Pheno[[1]]
  
  RSquared<-rep(NA,ncol(Matrix)) 
  Beta<-rep(NA,ncol(Matrix))
  P<-rep(NA,ncol(Matrix)) 
  for (i in 1:ncol(Matrix)) {
    mod<-lm(Subset[,i]~Pheno)
    RSquared[i]<-summary(mod)$r.squared
    Beta[i]<-summary(mod)$coefficients[2,1]
    P[i]<-summary(mod)$coefficients[2,4]
  }
  
  return(data.table(Metabolite=colnames(Matrix),Beta,RSquared,P,FDR=p.adjust(P,method="fdr"),
                    Bonferroni=p.adjust(P,method="bonferroni")))
}


MeanLS_Day04<-LM_Met_LS(Matrix,"MeanLS","Day04")
MeanLS_Day04[MeanLS_Day04$FDR<0.2,]
MeanLS_Day10<-LM_Met_LS(Matrix,"MeanLS","Day10")
MeanLS_Day10[MeanLS_Day10$FDR<0.2,]
MeanLS_Day24<-LM_Met_LS(Matrix,"MeanLS","Day24")
MeanLS_Day24[MeanLS_Day24$FDR<0.2,]
MeanLS_Day45<-LM_Met_LS(Matrix,"MeanLS","Day45")
MeanLS_Day45[MeanLS_Day45$FDR<0.2,]

logAlpha_Day04<-LM_Met_LS(Matrix,"logAlpha","Day04")
logAlpha_Day04[logAlpha_Day04$FDR<0.2,]
logAlpha_Day10<-LM_Met_LS(Matrix,"logAlpha","Day10")
logAlpha_Day10[logAlpha_Day10$FDR<0.2,]
logAlpha_Day24<-LM_Met_LS(Matrix,"logAlpha","Day24")
logAlpha_Day24[logAlpha_Day24$FDR<0.2,]
logAlpha_Day45<-LM_Met_LS(Matrix,"logAlpha","Day45")
logAlpha_Day45[logAlpha_Day45$FDR<0.2,]

Beta_Day04<-LM_Met_LS(Matrix,"Beta","Day04")
Beta_Day04[Beta_Day04$FDR<0.2,]
Beta_Day10<-LM_Met_LS(Matrix,"Beta","Day10")
Beta_Day10[Beta_Day10$FDR<0.2,]
Beta_Day24<-LM_Met_LS(Matrix,"Beta","Day24")
Beta_Day24[Beta_Day24$FDR<0.2,]
Beta_Day45<-LM_Met_LS(Matrix,"Beta","Day45")
Beta_Day45[Beta_Day45$FDR<0.2,]

######

### Correlation matrices of metabolite levels vs demography

## Collapse biological replicates 
Agg<-aggregate(Matrix,by=list(Genotype=SampleInfo$Genotype,Age=SampleInfo$Age),FUN=mean)

Agg4<-filter(Agg,Age=="Day04")

Agg4Matrix<-select(Agg4,-(Genotype:Age))
Agg4Matrix<-as.matrix(Agg4Matrix)
rownames(Agg4Matrix)<-Agg4$Genotype

Cor4<-corbetw2mat(PhenoMatrix,Agg4Matrix,what="all")

Agg10<-filter(Agg,Age=="Day10")

Agg10Matrix<-select(Agg10,-(Genotype:Age))
Agg10Matrix<-as.matrix(Agg10Matrix)
rownames(Agg10Matrix)<-Agg10$Genotype

Cor10<-corbetw2mat(PhenoMatrix,Agg10Matrix,what="all")

Agg24<-filter(Agg,Age=="Day24")

Agg24Matrix<-select(Agg24,-(Genotype:Age))
Agg24Matrix<-as.matrix(Agg24Matrix)
rownames(Agg24Matrix)<-Agg24$Genotype

Cor24<-corbetw2mat(PhenoMatrix,Agg24Matrix,what="all")

Agg45<-filter(Agg,Age=="Day45")
Pheno45<-filter(Pheno,Genotype!="Ral_796")

Pheno45Matrix<-select(Pheno45,-(Genotype:Block))
Pheno45Matrix<-as.matrix(Pheno45Matrix)
rownames(Pheno45Matrix)<-Pheno45$Genotype

Agg45Matrix<-select(Agg45,-(Genotype:Age))
Agg45Matrix<-as.matrix(Agg45Matrix)
rownames(Agg45Matrix)<-Agg45$Genotype

Cor45<-corbetw2mat(Pheno45Matrix,Agg45Matrix,what="all")

##########################################################

# Age-trajectory

load("CenterScaledData.RData")
Matrix<-CenterScaleBySample

SampleInfo<-fread("Metabolomics1_Randomized.txt")
SampleInfo<-SampleInfo[SampleID!="Sample181",]
SampleInfo<-select(SampleInfo,ExpBlock,Genotype,MetBatch,Age,SampleID)
SampleInfo[Age=="Day04", AgeNum:=4]
SampleInfo[Age=="Day10", AgeNum:=10]
SampleInfo[Age=="Day24", AgeNum:=24]
SampleInfo[Age=="Day45", AgeNum:=45]
SampleInfo[Age=="Day69", AgeNum:=69]
SampleInfo[Age=="Day80", AgeNum:=80]

Genotypes<-sort(unique(SampleInfo$Genotype))

## Remove samples from day 69 and 80 before fitting linear models

Index<-SampleInfo$Age!="Day69" & SampleInfo$Age!="Day80"
Matrix<-Matrix[Index,]
SampleInfo<-SampleInfo[Index,]

############

Wide<-data.table(Matrix)
Wide<-data.table(SampleInfo,Wide)

Beta_LM<-matrix(NA,nrow=length(Genotypes),ncol=ncol(Matrix))
rownames(Beta_LM)<-Genotypes
colnames(Beta_LM)<-colnames(Matrix)

for (i in 1:ncol(Matrix)) {
  for (j in 1:length(Genotypes)) {
    Data<-filter(Wide,Genotype==Genotypes[j])
    Time<-Data$AgeNum
    LM<-lm(Data[[i+6]]~Time)
    
    Beta_LM[j,i]<-summary(LM)$coefficients[2,1]
  }
}

Beta_LM
write.table(Beta_LM,"AgeTrajectory.txt",quote=F,sep="\t")

Shapiro<-apply(Beta_LM,2,function(x) shapiro.test(x)$p.value)
Shapiro[Shapiro<0.05]
length(Shapiro[Shapiro<0.05])

############

Pheno<-fread("Pheno20.txt")
Pheno<-select(Pheno,Genotype,Block,MeanLS,LogMAlpha,MBeta)
colnames(Pheno)[4:5]<-c("logAlpha","Beta")

PhenoMatrix<-select(Pheno,-(Genotype:Block))
PhenoMatrix<-as.matrix(PhenoMatrix)
rownames(PhenoMatrix)<-Pheno$Genotype

## Linear model of metabolite trajectories vs demograpy

LS<-Pheno$MeanLS
RSquared<-rep(NA,ncol(Matrix))
Beta<-rep(NA,ncol(Matrix))
P<-rep(NA,ncol(Matrix))

for (i in 1:ncol(Matrix)) {
  mod<-lm(Beta_LM[,i]~LS)
  RSquared[i]<-summary(mod)$r.squared
  Beta[i]<-summary(mod)$coefficients[2,1]
  P[i]<-summary(mod)$coefficients[2,4]
}

AgeTrajectory_LS<-data.table(Metabolite=colnames(Matrix),Beta,RSquared,P,FDR=p.adjust(P,method="fdr"))
filter(AgeTrajectory_LS,FDR<0.2)

#####

logAlpha<-Pheno$logAlpha
RSquared<-rep(NA,ncol(Matrix))
Beta<-rep(NA,ncol(Matrix))
P<-rep(NA,ncol(Matrix))

for (i in 1:ncol(Matrix)) {
  mod<-lm(Beta_LM[,i]~logAlpha)
  RSquared[i]<-summary(mod)$r.squared
  Beta[i]<-summary(mod)$coefficients[2,1]
  P[i]<-summary(mod)$coefficients[2,4]
}

AgeTrajectory_logAlpha<-data.table(Metabolite=colnames(Matrix),Beta,RSquared,P,FDR=p.adjust(P,method="fdr"))
filter(AgeTrajectory_logAlpha,FDR<0.2)

#####

RateOfAging<-Pheno$Beta
RSquared<-rep(NA,ncol(Matrix))
Beta<-rep(NA,ncol(Matrix))
P<-rep(NA,ncol(Matrix))

for (i in 1:ncol(Matrix)) {
  mod<-lm(Beta_LM[,i]~RateOfAging)
  RSquared[i]<-summary(mod)$r.squared
  Beta[i]<-summary(mod)$coefficients[2,1]
  P[i]<-summary(mod)$coefficients[2,4]
}

AgeTrajectory_Beta<-data.table(Metabolite=colnames(Matrix),Beta,RSquared,P,FDR=p.adjust(P,method="fdr"))
filter(AgeTrajectory_Beta,FDR<0.2)

CorTraj<-corbetw2mat(PhenoMatrix,Beta_LM,what="all")

################

## Permutation testing if the correlation I see is unlikely by chance

Permutation<-function(Vector) {
  Num<-100 # Number of permutation
  
  P<-matrix(NA,nrow=Num,ncol=ncol(Matrix))
  colnames(P)<-colnames(Matrix)
  FDR<-matrix(NA,nrow=Num,ncol=ncol(Matrix))
  colnames(FDR)<-colnames(Matrix)
  
  for (i in 1:Num) {
    Fake<-sample(Vector)
    for (j in 1:ncol(Matrix)) {
      mod<-lm(Beta_LM[,j]~Fake)
      P[i,j]<-summary(mod)$coefficients[2,4]
      FDR[i,]=p.adjust(P[i,],method="fdr")
    }
  }
  
  FDR0.2<-apply(FDR,1,function(x) sum(x<0.2))
  plot(density(FDR0.2))
  abline(v=5,col="red")
  
  MoreExtreme<-sum(FDR0.2>=5)/Num
  return(MoreExtreme) 
}

Permutation(LS)
Permutation(logAlpha)
Permutation(RateOfAging)

#####################################################################

# Visualization

CorMeanLS<-rbind(Cor4[1,],Cor10[1,],Cor24[1,],Cor45[1,],CorTraj[1,])
rownames(CorMeanLS)<-c("Day04","Day10","Day24","Day45","AgeTrajectory")

CorAlpha<-rbind(Cor4[2,],Cor10[2,],Cor24[2,],Cor45[2,],CorTraj[2,])
rownames(CorAlpha)<-c("Day04","Day10","Day24","Day45","AgeTrajectory")

CorBeta<-rbind(Cor4[3,],Cor10[3,],Cor24[3,],Cor45[3,],CorTraj[3,])
rownames(CorBeta)<-c("Day04","Day10","Day24","Day45","AgeTrajectory")


h1<-Heatmap(CorMeanLS,name="Correlation Coefficient",
            show_row_names=TRUE,show_column_names=FALSE,
            show_row_dend=FALSE,show_column_dend=FALSE,
            cluster_rows = FALSE,cluster_columns = TRUE)
h2<-Heatmap(CorAlpha,name="Correlation Coefficient",
            show_row_names=TRUE,show_column_names=FALSE,
            show_row_dend=FALSE,show_column_dend=FALSE,
            cluster_rows = FALSE,cluster_columns = TRUE)
h3<-Heatmap(CorBeta,name="Correlation Coefficient",
            show_row_names=TRUE,show_column_names=FALSE,
            show_row_dend=FALSE,show_column_dend=FALSE,
            cluster_rows = FALSE,cluster_columns = TRUE)

pdf("Met_MeanLS.pdf", width=10, height=2)
h1
dev.off()

pdf("Met_Alpha.pdf", width=10, height=2)
h2
dev.off()

pdf("Met_Beta.pdf", width=10, height=2)
h3
dev.off()




