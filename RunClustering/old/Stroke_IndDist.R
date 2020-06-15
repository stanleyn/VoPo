library('flowCore')
library('foreach')
library('doParallel')
library('iterators')
library('Rtsne')

source('Setabr/runRepMetaclust.R')
source('Setabr/runClassif.R')
source('Setabr/vizAtlas.R')

###############################
##Step 1: Put in filename, info, etc.
##################################

#Modify here!!! Put the path to your FCS files
source('~/Stroke/Process_Stroke_June.R')
FileNames=FNames
###############################################

MN=readRDS('Processed/MN_Stroke')
ToUse=readRDS('Processed/ToUse_Stroke')
Meta_Stroke=readRDS('Processed/Meta_Stroke')

###############################
#Clustering
###############################

#AUCBoot=c()
#for(s in 1:100){
#Build=runRepMetaclust(50,50,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)

#print('done with clustering')

####################
#Classification
################
#IterNumClus=Build$IterNumClus
#FuncDF=Build[[2]]
#colnames(FuncDF)=1:ncol(FuncDF)

#AUCs=runClassif(FuncDF,factor(Meta_Stroke$Class),40,IterNumClus,0.2,1,as.character(Meta_Stroke$Subject),35)

#print('done with classification')
#AUCBoot=c(AUCBoot,AUCs)
#}
###################################################
AUCBase=c()
for(s in 1:100){
Build=runRepMetaclust(1,50,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)

print('done with clustering')

####################
#Classification
################
IterNumClus=Build$IterNumClus
FuncDF=Build[[2]]
colnames(FuncDF)=1:ncol(FuncDF)

AUCs=runClassif(FuncDF,factor(Meta_Stroke$Class),IterNumClus[1],IterNumClus,0.2,1,as.character(Meta_Stroke$Subject),35)

print('done with classification')
AUCBase=c(AUCBase,AUCs)
}

