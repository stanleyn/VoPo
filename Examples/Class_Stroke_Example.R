#purpose: generate a distribution of AUC
library('ggplot2')
library('reshape2')
source('VoPo_main/runClassif.R')
source('VoPo_main/getFrequencyFeature.R')

#get processed VoPo result
Build_Stroke=readRDS('Processed/Build_Stroke.rds')

#load metadata about samples, etc.
Meta_Stroke=readRDS('Processed/Meta_Stroke.rds')

#use your vector of filenames for the rows of your feature matrix
NameVec=Meta_Stroke$FileNames

#extract features from VoPo
FuncDF=getFrequencyFeature(Build_Stroke,NameVec)

ClAcc=runClassif(FuncDF,factor(Meta_Stroke$Class),40,Build_Stroke$IterNumClus,0.5,100,as.character(Meta_Stroke$Subject),35)

