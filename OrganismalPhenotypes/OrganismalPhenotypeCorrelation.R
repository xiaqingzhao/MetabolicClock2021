library(data.table)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(corrplot)

Pheno<-fread("Pheno_20.txt")

Matrix<-as.matrix(Pheno[,c(3,5,6,14,15,16,17,18,19,20,21,25)])
rownames(Matrix)<-Pheno$Genotype
colnames(Matrix)<-c("MeanLS",":log(alpha)",":beta","RepOutput_Day08","RepOutput_Day12",
                    "Activity_Week1","Activity_Week2","Activity_Week3","Activity_Week4","Activity_Week5","Activity_Week6","ActivityDecline")

Cor<-cor(Matrix,method="spearman",use="pairwise.complete.obs")

pdf("PhenotypeCor_Genotype20.pdf", width=5, height=5)
corrplot(Cor, method = "circle",type="lower",diag=FALSE)
dev.off()

##########

shapiro.test(Pheno$MeanLS)
shapiro.test(Pheno$LogMAlpha)
shapiro.test(Pheno$MBeta)
shapiro.test(Pheno$RepOutput_Day08)
shapiro.test(Pheno$RepOutput_Day12)
shapiro.test(Pheno$Week1)
shapiro.test(Pheno$Week2)
shapiro.test(Pheno$Week3)
shapiro.test(Pheno$Week4)
shapiro.test(Pheno$Week5)
shapiro.test(Pheno$Week6)
shapiro.test(Pheno$Slope)

##########

cor.test(Pheno$MeanLS,Pheno$RepOutput_Day08, method="spearman")
cor.test(Pheno$MeanLS,Pheno$RepOutput_Day12, method="spearman")
cor.test(Pheno$MeanLS,Pheno$Week1, method="spearman")
cor.test(Pheno$MeanLS,Pheno$Week2, method="spearman")
cor.test(Pheno$MeanLS,Pheno$Week3, method="spearman")
cor.test(Pheno$MeanLS,Pheno$Week4, method="spearman")
cor.test(Pheno$MeanLS,Pheno$Week5, method="spearman")
cor.test(Pheno$MeanLS,Pheno$Week6, method="spearman")
cor.test(Pheno$MeanLS,Pheno$Slope, method="spearman")

cor.test(Pheno$LogMAlpha,Pheno$RepOutput_Day08, method="spearman")
cor.test(Pheno$LogMAlpha,Pheno$RepOutput_Day12, method="spearman")
cor.test(Pheno$LogMAlpha,Pheno$Week1, method="spearman")
cor.test(Pheno$LogMAlpha,Pheno$Week2, method="spearman")
cor.test(Pheno$LogMAlpha,Pheno$Week3, method="spearman")
cor.test(Pheno$LogMAlpha,Pheno$Week5, method="spearman")
cor.test(Pheno$LogMAlpha,Pheno$Slope, method="spearman")

cor.test(Pheno$MBeta,Pheno$RepOutput_Day08, method="spearman")
cor.test(Pheno$MBeta,Pheno$RepOutput_Day12, method="spearman")
cor.test(Pheno$MBeta,Pheno$Week1, method="spearman")
cor.test(Pheno$MBeta,Pheno$Week2, method="spearman")
cor.test(Pheno$MBeta,Pheno$Week3, method="spearman")
cor.test(Pheno$MBeta,Pheno$Week4, method="spearman")
cor.test(Pheno$MBeta,Pheno$Week5, method="spearman")
cor.test(Pheno$MBeta,Pheno$Week6, method="spearman")
cor.test(Pheno$MBeta,Pheno$Slope, method="spearman")

