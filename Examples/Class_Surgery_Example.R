#purpose: generate a distribution of AUC
source('VoPo_main/runClassif.R')
library('reshape2')
library('ggplot2')
source('VoPo_main/getFrequencyFeature.R')

#processed VoPo result
Build_Surgery=readRDS('Processed/Build_Surgery.rds')

#load metadata about the samples, etc
Meta_Surgery=readRDS('Processed/Meta_Surgery.rds')

#use the vector of file names as the rownames for your feature matrix
NameVec=Meta_Surgery$FileNames

#Extract the features from the VoPo object
FuncDF=getFrequencyFeature(Build_Surgery,NameVec)

#get repeated distribution
ClAcc=runClassif(FuncDF,Meta_Surgery$Class,40,Build_Surgery$IterNumClus,0.5,100,as.character(Meta_Surgery$Subject),35)


