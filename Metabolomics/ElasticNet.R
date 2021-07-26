library(data.table)
library(tidyverse)
library(caret)
library(gridExtra)

load("CenterScaledData.RData")
Matrix<-CenterScaleBySample
ScaledMatrix<-scale(Matrix)

SampleInfo<-fread("Metabolomics1_Randomized.txt")
SampleInfo<-SampleInfo[SampleID!="Sample181",]
SampleInfo<-select(SampleInfo,ExpBlock,Genotype,MetBatch,Age,SampleID)
SampleInfo[Age=="Day04", AgeNum:=4]
SampleInfo[Age=="Day10", AgeNum:=10]
SampleInfo[Age=="Day24", AgeNum:=24]
SampleInfo[Age=="Day45", AgeNum:=45]
SampleInfo[Age=="Day69", AgeNum:=69]
SampleInfo[Age=="Day80", AgeNum:=80]

SampleAge<-SampleInfo$AgeNum

#########################################

## 1. Fix a value of alpha (turned out to be 0)

EN<-function(Omics,Pheno) {
  
  ## Split data into training and test sets
  set.seed(1000)
  index<-sample(1:nrow(Omics),round(nrow(Omics)*0.8))
  Omics_train<-Omics[index,]
  Omics_test<-Omics[-index,]
  Pheno_train<-Pheno[index]
  Pheno_test<-Pheno[-index]
  
  ## Tune grid
  tuneGrid = expand.grid(
    alpha=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
    # alpha=0,
    lambda=c(seq(0.000001,1,length=100),seq(1,200,length=201))
  )
  
  ## Model fitting 
  model<-train(
    Omics_train, 
    Pheno_train,
    tuneGrid = tuneGrid,
    method = "glmnet", 
    trControl = trainControl(
      method = "cv", 
      number = 5,  ## 5-fold cross validation
      verboseIter = TRUE
    )
  )
  
  error<-Pheno_test-predict(model,newdata=Omics_test)
  RMSE<-sqrt(mean(error^2)) 
  
  return(list(Model=model,test_RMSE=RMSE))
}

InitialEN<-EN(Matrix,SampleAge)
InitialEN
pdf("ElasticNetParameter.pdf", width=10, height=4)
plot(InitialEN[[1]])
dev.off()

## 2. On the full data, use cross validation to learn a value of lambda (turned out to be 5)

EN.full<-function(Omics,Pheno,tuneGrid) {
  model <- train(
    Omics, 
    Pheno,
    tuneGrid = tuneGrid,
    method = "glmnet",
    trControl = trainControl(
      method = "cv", 
      number = 5, 
      verboseIter = TRUE
    )
  )
}

OptLambda<-rep(NA,100)
for (i in 1:100) {
  set.seed(i)
  tuneGrid = expand.grid(alpha=0,lambda=c(seq(0.000001,1,length=100),seq(1,1000,length=1001)))
  FullEN<-EN.full(Matrix,SampleAge,tuneGrid)
  OptLambda[i]<-FullEN$bestTune$lambda
}

tuneGrid = expand.grid(alpha=0,lambda=c(seq(0.000001,1,length=100),seq(1,1000,length=1001)))
FullEN<-EN.full(Matrix,SampleAge,tuneGrid)
FullEN$bestTune$lambda

## 3. Get prediction on each sample by: 
#### Removing one sample from the dataset
#### Fitting a model on the remaining samples, using the fixed alpha and lambda values from 1 and 2
#### Getting a prediction on the left-out sample

EN_LOO<-function(Omics,Pheno) {
  
  PredAge<-rep(NA,nrow(Omics))
  
  for (i in 1:nrow(Omics)) {
    index<-i
    Omics_train<-Omics[-c(index),]
    Omics_test<-Omics[index,,drop=FALSE]
    Pheno_train<-Pheno[-c(index)]
    Pheno_test<-Pheno[index]
    
    tuneGrid=expand.grid(alpha=0,lambda=5)
    
    model<-train(
      Omics_train, 
      Pheno_train,
      tuneGrid = tuneGrid,
      method = "glmnet", 
      trControl = trainControl(
        method = "cv", 
        number = 5,  ## 5-fold cross validation
        verboseIter = TRUE
      )
    )
    
    PredAge[i]<-predict(model,newdata=Omics_test)
    
  }
  
  return(PredAge)
}

EN_LOO_Result<-EN_LOO(Matrix,SampleAge)
EN_LOO_Result<-data.frame(SampleInfo,PredAge=EN_LOO_Result)
write.table(EN_LOO_Result,"ElasticNetPrediction.txt",quote=F,row.names=F,sep="\t")

################

## Residuals

library(ggforce)

Residuals<-data.table(Genotype=SampleInfo$Genotype,Age=SampleInfo$Age,AgeNum=SampleInfo$AgeNum,
                      Residual=(EN_LOO_Result$PredAge-EN_LOO_Result$AgeNum))

Pheno<-fread("Pheno20.txt")
Pheno<-select(Pheno,Genotype,Block,MeanLS,LogMAlpha,MBeta,LogMiu45,LogMiu60)
colnames(Pheno)[4:7]<-c("logAlpha","Beta","logMiu45","logMiu60")
PhenoMatrix<-select(Pheno,-(Genotype:Block))
PhenoMatrix<-as.matrix(PhenoMatrix)
rownames(PhenoMatrix)<-Pheno$Genotype 

DT<-merge(Residuals,Pheno,by="Genotype")
LongDT<-gather(DT,Phenotype,Measurement,MeanLS:logMiu60)
LongDT$Phenotype<-factor(LongDT$Phenotype,levels=c("MeanLS","logAlpha","Beta","logMiu45","logMiu60"))

LongDT_NoLateAges<-filter(LongDT,Age!="Day69" & Age!="Day80")

# Plot<-ggplot(na.omit(LongDT_NoLateAges),aes(x=Measurement,y=Residual))+geom_point()+
#   theme_bw()+geom_smooth(method="lm",col="red")+facet_grid(Age~Phenotype,scales='free_x')+
#   ylab("Age Acceleration")+theme(axis.text.x=element_text(size=6),axis.title.x=element_blank())

LongDT_NoLateAges$Phenotype2<-factor(LongDT_NoLateAges$Phenotype,
                                     labels=c("MeanLS","'log('*alpha*')'","beta","'log('*mu*')'[45]","'log('*mu*')'[60]"))

Plot<-ggplot(na.omit(LongDT_NoLateAges),aes(x=Measurement,y=Residual))+geom_point()+
  theme_bw()+geom_smooth(method="lm",col="red")+
  facet_grid(Age~Phenotype2,scales='free_x',labeller=label_parsed)+
  ylab("Age Acceleration")+theme(axis.text.x=element_text(size=6),axis.title.x=element_blank())

pdf("AgeAcceleration.pdf",width=8,height=6)
Plot
dev.off()

## Correlation

Day45AA<-DT[DT$Age=="Day45",]
cor.test(Day45AA$Residual,Day45AA$MeanLS,method="spearman")
cor.test(Day45AA$Residual,Day45AA$logAlpha,method="spearman")
cor.test(Day45AA$Residual,Day45AA$logMiu45,method="spearman")
cor.test(Day45AA$Residual,Day45AA$logMiu60,method="spearman")


#########################

## Plot with true training and testing

EN.truetest<-function(Omics,Pheno,tuneGrid) {
  
  ## Split data into training and test sets
  set.seed(1001)
  index<-sample(1:nrow(Omics),round(nrow(Omics)*0.8))
  Omics_train<-Omics[index,]
  Omics_test<-Omics[-index,]
  Pheno_train<-Pheno[index]
  Pheno_test<-Pheno[-index]
  
  ## Model fitting 
  set.seed(4321)
  model<-train(
    Omics_train, 
    Pheno_train,
    tuneGrid = tuneGrid,
    method = "glmnet", 
    trControl = trainControl(
      method = "cv", 
      number = 5,  
      verboseIter = TRUE
    )
  )
  
  Pred_test<-predict(model,newdata=Omics_test)
  error<-Pheno_test-Pred_test
  RMSE<-sqrt(mean(error^2)) 
  
  Pred_train<-predict(model,newdata=Omics_train)
  
  return(list(Model=model,test_RMSE=RMSE,TrueVSPred_Test=cbind(Pheno_test,Pred_test),
              TrueVSPred_Train=cbind(Pheno_train,Pred_train)))
}

tuneGrid = expand.grid(alpha=0,lambda=4)

ElasticNet<-EN.truetest(Matrix,SampleAge,tuneGrid)

e<-ggplot(data.frame(ElasticNet[[3]]),aes(x=Pheno_test,y=Pred_test))+geom_point(size=0.5)+geom_abline(slope=1,intercept=0,col="red",linetype = "dashed")+
  scale_x_continuous(limits = c(0, 80))+scale_y_continuous(limits = c(-4, 80))+theme_bw()+
  xlab("Real Sample Age")+ylab("Predicted Sample Age")+geom_smooth(method="lm",col="blue")+ggtitle("Test set")
f<-ggplot(data.frame(ElasticNet[[4]]),aes(x=Pheno_train,y=Pred_train))+geom_point(size=0.5)+geom_abline(slope=1,intercept=0,col="red",linetype = "dashed")+
  scale_x_continuous(limits = c(0, 80))+scale_y_continuous(limits = c(-4, 80))+theme_bw()+
  xlab("Real Sample Age")+ylab("Predicted Sample Age")+geom_smooth(method="lm",col="blue")+ggtitle("Training set")

pdf("RidgeRegression_TrainingTesting.pdf",width=8,height=4)
grid.arrange(f,e,ncol=2)
dev.off()

cor.test(ElasticNet[[3]][,1],ElasticNet[[3]][,2])

#################


