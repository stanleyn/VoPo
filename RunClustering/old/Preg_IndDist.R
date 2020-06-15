#Purpose: This is the steps needed to run bootsrapped BL1 on Kehlet
#Date Created: July 25, 2018
######################################################

#dependencies
library('flowCore')
library('foreach')
library('doParallel')
library('iterators')
library('Rtsne')

source('Setabr/runRepMetaclust.R')
source('Setabr/runClassif.R')
source('Setabr/vizAtlas.R')

########################################
##Step 1: Put in filenames, info, etc.
########################################

#Modify here!!! Put the path to your FCS files
source('~/Pregnancy/Tools/Process_PP.R')
################################################

MN=readRDS('Processed/MN_Pregnancy')
ToUse=readRDS('Processed/ToUse_Pregnancy')
Meta_Preg=readRDS('Processed/Meta_Pregnancy')

#AUCSingle=c()

#for(s in 1:100){

#################
#clustering
################

#Build=runRepMetaclust(1,30,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)

#print('clustering done')

################
#Classification
###############

#print('running 30 trials of the classification pipeline')

#FuncDF=Build[[2]]
#IterNumClus=Build$IterNumClus
#ForC=1:ncol(FuncDF)
#NewName=paste(colnames(FuncDF),ForC,sep='_')
#colnames(FuncDF)=NewName
#AUCs=runClassif(FuncDF,Meta_Preg$Class,IterNumClus[1],IterNumClus,0.5,1,as.character(Meta_Preg$Subject),35)

#AUCSingle=c(AUCSingle,AUCs)
#print('done with classification')
#}

AUCAll=c()

for(s in 1:100){

#################
#clustering
################

Build=runRepMetaclust(50,30,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)

print('clustering done')

################
#Classification
###############

print('running 30 trials of the classification pipeline')

FuncDF=Build[[2]]
IterNumClus=Build$IterNumClus
ForC=1:ncol(FuncDF)
NewName=paste(colnames(FuncDF),ForC,sep='_')
colnames(FuncDF)=NewName
AUCs=runClassif(FuncDF,Meta_Preg$Class,10,IterNumClus,0.5,1,as.character(Meta_Preg$Subject),35)

AUCAll=c(AUCAll,AUCs)
print('done with classification')
}

