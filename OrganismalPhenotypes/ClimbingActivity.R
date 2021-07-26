library(data.table)
library(tidyverse)
library(lme4)
library(ggforce)

Data<-fread("ClimbingActivity20_Raw.txt")

## Collapse flies within the same vial 
Climb_MeanAbs<-Data[,.(mean(AdjFlyHeight)),by=.(Block,Genotype,Week,Chamber)]
Climb_SEAbs<-Data[,.(sd(AdjFlyHeight)/sqrt(.N)),by=.(Block,Genotype,Week,Chamber)]

Climb_MeanP<-Data[,.(mean(AdjFlyHeight/AdjMaxHeight)),by=.(Block,Genotype,Week,Chamber)]
Climb_SEP<-Data[,.(sd(AdjFlyHeight/AdjMaxHeight)/sqrt(.N)),by=.(Block,Genotype,Week,Chamber)]

Climb<-data.frame(Block=Climb_MeanAbs$Block,Genotype=Climb_MeanAbs$Genotype, 
                  Week=Climb_MeanAbs$Week,Mean=Climb_MeanAbs$V1,SE=Climb_SEAbs$V1)
ClimbP<-data.frame(Block=Climb_MeanAbs$Block,Genotype=Climb_MeanAbs$Genotype, 
                   Week=Climb_MeanAbs$Week,Mean=Climb_MeanP$V1,SE=Climb_SEP$V1)

#########

fit<-lm(Mean~Block+Genotype+as.character(Week),data=Climb)
anova(fit)
#plot(fit)

fit<-lm(Mean~Block+Genotype+as.character(Week),data=ClimbP)
anova(fit)
#plot(fit)

## Using the percentage climbed helps to reduce block effect. So proceed with percentages
## Model fitting seems good for both absolute height and percentage height

########

## Visualization 

Plot<-ggplot(ClimbP,aes(y=Mean,x=as.character(Week)))+geom_boxplot()+theme_bw()
Plot+facet_wrap('Genotype')

########

## Collapse biological replicates 
ClimbSummary<-ClimbP %>% group_by(Block,Genotype,Week) %>% summarise(PHeight=mean(Mean))

fit<-lm(PHeight~Block+Genotype+as.character(Week),data=ClimbSummary)
anova(fit)

########

## Turn summary data into wide format

Climb_Wide<-ClimbSummary
Climb_Wide$Age<-paste("Week",ClimbSummary$Week,sep="")
Climb_Wide<-select(Climb_Wide,-Week)
Climb_Wide<-spread(Climb_Wide,Age,PHeight)
Climb_Wide<-data.table(Climb_Wide)
Climb_Wide<-Climb_Wide[order(Genotype)]

########

## Baseline climbing ability and rate of climbing ability decrease 

Intercept<-rep(NA,nrow(Climb_Wide))
InterceptAt3<-rep(NA,nrow(Climb_Wide))
InterceptAt5<-rep(NA,nrow(Climb_Wide))
Slope<-rep(NA,nrow(Climb_Wide))
for (i in 1:nrow(Climb_Wide)) {
  if (is.na(Climb_Wide[i,Week5])==FALSE | is.na(Climb_Wide[i,Week6])==FALSE) {
    Week<-c(3,4,5,6)
    mod<-lm(Climb_Wide[i,c(Week3,Week4,Week5,Week6)]~Week)
    Intercept[i]<-mod$coefficients[1]
    Slope[i]<-mod$coefficients[2]
    InterceptAt3[i]<-predict(mod)[1]
    InterceptAt5[i]<-predict(mod)[3]
  }
  else {}
}
Summary_LM<-data.table(Climb_Wide[,1:2],Intercept,InterceptAt3,InterceptAt5,Slope)


Climb_Result<-data.table(Climb_Wide,Summary_LM[,3:6])

Pairs<-as.matrix(Climb_Result[,-c(1,2)])
row.names(Pairs)<-Climb_Wide$Genotype
pairs(Pairs,pch=20)

write.table(Climb_Result,"Activity_Genotype20.txt",quote=F,sep="\t",row.names=F)

########

## Plot of example strains

Exp<-ClimbP[ClimbP$Genotype %in% c("Ral_441","Ral_355","Ral_136"),]
Exp$Genotype=factor(Exp$Genotype, levels=c('Ral_441','Ral_355','Ral_136'))

Plot<-ggplot(Exp,aes(y=Mean,x=as.character(Week)))+geom_boxplot()+theme_bw()+xlab("Week")+ylab("Activity Level")
pdf("ActivityLevelExample.pdf", width=8, height=3)
Plot+facet_wrap('Genotype')
dev.off()

##############

## Visualize variation of climbing over age

Fill<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')
ClimbSummary<-arrange(ClimbSummary,Genotype,Week)
ClimbSummary<-filter(ClimbSummary, Genotype!="Ral_101" & Genotype!="Ral_38" & Genotype!="Ral_385")

Plot<-ggplot(ClimbSummary, aes(x=Week,y=PHeight,group=Genotype,colour=Genotype))+
  geom_point()+theme_bw()+xlab("Week")+ylab("Activity Level")+
  scale_color_manual(values=Fill)+theme(legend.position="none")

GenoFill<-data.frame(Genotype=ClimbSummary$Genotype,GenoCol=ggplot_build(Plot)$data[[1]]$colour)
GenoFill<-unique(GenoFill)

ClimbVar<-Climb_Result[,c(1,2,9,12)]
ClimbVar<-merge(ClimbVar,GenoFill)
ClimbVar<-ClimbVar[complete.cases(ClimbVar),]

for (i in 1:nrow(ClimbVar)) {
  Plot<-Plot+geom_abline(slope=ClimbVar$Slope[i], intercept=ClimbVar$Intercept[i], colour=ClimbVar$GenoCol[i])+
    scale_x_continuous(limits=c(3, 6))
}

pdf("VariationClimbing.pdf",width=4,height=4)
Plot
dev.off()



