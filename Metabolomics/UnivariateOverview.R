library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(car)
library(lme4)

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
Age.ordered<-factor(SampleInfo$Age,levels=c("Day04","Day10","Day24","Day45","Day69","Day80"),ordered=TRUE)


###############

## Remove samples from day 80 before fitting linear models

Index<-SampleInfo$Age!="Day80"
Matrix<-Matrix[Index,]
SampleInfo<-SampleInfo[Index,]
Age.ordered<-Age.ordered[Index]

###############

Wide<-data.table(Matrix)
Wide<-data.table(SampleInfo,Wide)
Long<-gather(Wide,Compound,Measurement,-(ExpBlock:AgeNum))

#######################

## Linear Mixed Model 

GenotypeVar<-rep(NA,ncol(Matrix))
MetBatchVar<-rep(NA,ncol(Matrix))
ExpBlockVar<-rep(NA,ncol(Matrix))
ResidualVar<-rep(NA,ncol(Matrix))
AgeVar<-rep(NA,ncol(Matrix))
AgeP<-rep(NA,ncol(Matrix))

for (i in 1:ncol(Matrix)) {
  DF<-data.frame(Measurement=Wide[[i+6]],Age=Wide$AgeNum,Genotype=Wide$Genotype,ExpBlock=Wide$ExpBlock,MetBatch=Wide$MetBatch)
  FM<-lmer(Measurement~Age+(1|Genotype)+(1|ExpBlock)+(1|MetBatch),REML=FALSE,data=DF)
  
  Var<-as.data.frame(VarCorr(FM))
  GenotypeVar[i]<-Var[1,4]
  MetBatchVar[i]<-Var[2,4]
  ExpBlockVar[i]<-Var[3,4]
  ResidualVar[i]<-Var[4,4]
  
  AgeVar[i]<-(cor(DF$Age,DF$Measurement))^2
  AgeP[i]<-Anova(FM)$'Pr(>Chisq)'
}

FM_Result<-data.table(Compound=colnames(Matrix),GenotypeVar=GenotypeVar,MetBatchVar=MetBatchVar,
                      ExpBlockVar=ExpBlockVar,ResidualVar=ResidualVar,AgeVar=AgeVar,AgeP=AgeP)
setDT(FM_Result)[, TotalVar:=GenotypeVar+MetBatchVar+ExpBlockVar+ResidualVar+AgeVar]
VarComponent<-mutate(FM_Result,Genotype=100*GenotypeVar/TotalVar,MetBatch=100*MetBatchVar/TotalVar,
                     ExpBlock=100*ExpBlockVar/TotalVar,Residual=100*ResidualVar/TotalVar,Age=100*AgeVar/TotalVar)

#### 

### Visualize variance partition 

FM_Result_Ordered<-FM_Result[order(TotalVar)]
VarPar<-select(FM_Result_Ordered,Compound,GenotypeVar,MetBatchVar,ExpBlockVar,ResidualVar,AgeVar)
LongVarPar<-gather(VarPar,Source,Variance,GenotypeVar:AgeVar)
LongVarPar<-LongVarPar %>% mutate(Compound=factor(Compound,levels=FM_Result_Ordered$Compound))
LongVarPar$Source<-factor(LongVarPar$Source, levels=c("AgeVar","GenotypeVar","ExpBlockVar","MetBatchVar","ResidualVar"))

pdf("VariancePartition.pdf", width=12, height=5)
ggplot(LongVarPar,aes(fill=Source,y=Variance,x=Compound))+
  geom_bar(position="stack",stat="identity")+theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x = element_text(angle=90))+xlab("Metabolite")+
  scale_fill_manual(values=c("#F8766D", "#00C19C", "#00B0F6","#BB9D00","#C77CFF"),labels=c("Age", "Genotype", "Experiment Block","Metabolomics Batch","Residual"))
dev.off()

#######################

## Linear model, with interaction between genotype and age

Age_Beta<-rep(NA,ncol(Matrix))
Age_P<-rep(NA,ncol(Matrix))
Genotype_P<-rep(NA,ncol(Matrix))
Interaction_P<-rep(NA,ncol(Matrix))

for (i in 1:ncol(Matrix)) {
  LM<-lm(Wide[[i+6]]~Wide$AgeNum*Wide$Genotype)
  Age_Beta[i]<-coef(LM)[2]
  Age_P[i]<-anova(LM)$'Pr(>F)'[1]
  Genotype_P[i]<-anova(LM)$'Pr(>F)'[2]
  Interaction_P[i]<-anova(LM)$'Pr(>F)'[3]
}

LM_Inter_Result<-data.table(Compound=colnames(Matrix),Age_Beta,
                            Age_P,Age_FDR=p.adjust(Age_P,method="fdr"),
                            Genotype_P,Genotype_FDR=p.adjust(Genotype_P,method="fdr"),
                            Interaction_P,Interaction_FDR=p.adjust(Interaction_P,method="fdr"))
LM_Inter_Result[Age_FDR<0.01] 
LM_Inter_Result[Age_FDR<0.01 & Age_Beta>0]
LM_Inter_Result[Age_FDR<0.01 & Age_Beta<0]
LM_Inter_Result[Genotype_FDR<0.01]
LM_Inter_Result[Interaction_FDR<0.01] ## Genotype-Age interaction seems to be present in a large number of metabolites

###################

## Venn diagram

All<-LM_Inter_Result$Compound
Age<-LM_Inter_Result[Age_FDR<0.01]$Compound
IncreaseWithAge<-LM_Inter_Result[Age_FDR<0.01 & Age_Beta>0]$Compound
DecreaseWithAge<-LM_Inter_Result[Age_FDR<0.01 & Age_Beta<0]$Compound
Genotype<-LM_Inter_Result[Genotype_FDR<0.01]$Compound
AgeGenotypeInteraction<-LM_Inter_Result[Interaction_FDR<0.01]$Compound

x<-list(
  Age=Age,
  Genotype=Genotype,Interaction=AgeGenotypeInteraction
)

y<-list(
  IncreaseWithAge=IncreaseWithAge,DecreaseWithAge=DecreaseWithAge,
  Genotype=Genotype, Interaction=AgeGenotypeInteraction
)

z<-list(All=All, Age=Age,
        Genotype=Genotype,Interaction=AgeGenotypeInteraction)

library(eulerr)

pdf("VennDiagram1.pdf",width=4,height=4)
plot(euler(x,shape="ellipse"),quantities=TRUE)
dev.off()

# plot(euler(z,shape="ellipse"),quantities=TRUE)
# plot(euler(z),quantities=TRUE)

# plot(venn(x))
# plot(venn(y))

# plot(euler(x,shape="ellipse"),quantities=TRUE)
# plot(euler(y,shape="ellipse"),quantities=TRUE)


