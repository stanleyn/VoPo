#purpose: generate a distribution of AUC
source('VoPo_main/runClassif.R')
library('ggplot2')
library('reshape2')
source('VoPo_main/getFrequencyFeature.R')


#processed VoPo result
Build_Preg=readRDS('Processed/Build_Pregnancy.rds')

#load metadata about samples, etc. 
Meta_Preg=readRDS('Processed/Meta_Pregnancy.rds')

#use your vector of filenames or a vector of names for the rows of your feature matrix
NameVec=Meta_Preg$FileNames

#extract features from VoPo object
FuncDF=getFrequencyFeature(Build_Preg,NameVec)

ClAcc=runClassif(FuncDF,Meta_Preg$Class,40,Build_Preg$IterNumClus,0.5,100,as.character(Meta_Preg$Subject),35)


