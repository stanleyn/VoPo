#Purpose: This is the steps needed to run bootsrapped BL2 on Kehlet
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
##Step 1: Put in filenames, info etc. #
########################################

#Modify here!!! Put the path to your FCS files
FileNames=list.files(path='~/FCS',pattern='.fcs',full.names=TRUE)
###################################

MN=readRDS('Processed/MN_Surgery')
ToUse=readRDS('Processed/ToUse_Surgery')
Meta_Surgery=readRDS('Processed/Meta_Surgery')

##Single Experiment##

#AUCSingle=c()
#for(i in 1:2){

############
#clustering
############
#Build=runRepMetaclust(1,30,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)

#print('clustering done!')

################
#classification
###############

#FuncDF=Build[[2]]
#IterNumClus=Build$IterNumClus
#make sure that each column of FuncDF has a different name
#ForC=1:ncol(FuncDF)
#NewName=paste(colnames(FuncDF),ForC,sep='_')
#colnames(FuncDF)=NewName

#print('running 30 trials of the classification pipeline')
#AUCs=runClassif(FuncDF,Meta_Surgery$Class,IterNumClus[1],IterNumClus,0.5,1,as.character(Meta_Surgery$Subject),35)
#AUCSingle=c(AUCSingle,AUCs)
#print('classification done!')
#}

AUCBoot=c()
for(i in 1:54){

############
#clustering
############
Build=runRepMetaclust(50,30,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)

print('clustering done!')

################
#classification
###############

FuncDF=Build[[2]]
IterNumClus=Build$IterNumClus
#make sure that each column of FuncDF has a different name
ForC=1:ncol(FuncDF)
NewName=paste(colnames(FuncDF),ForC,sep='_')
colnames(FuncDF)=NewName

print('running 30 trials of the classification pipeline')
AUCs=runClassif(FuncDF,Meta_Surgery$Class,10,IterNumClus,0.5,1,as.character(Meta_Surgery$Subject),35)
AUCBoot=c(AUCBoot,AUCs)
print('classification done!')
}

