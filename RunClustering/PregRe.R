#Purpose: This is the steps needed to run bootsrapped BL1 on Kehlet
#Date Created: July 25, 2018
######################################################

#dependencies
library('flowCore')
library('foreach')
library('doParallel')
library('iterators')
library('Rtsne')

source('VoPo_main/runRepMetaclust.R')
source('VoPo_main/runClassif.R')
source('VoPo_main/vizAtlas.R')

########################################
##Step 1: Put in filenames, info, etc.
########################################

#Modify here!!! Put the path to your FCS files
#FileNames=list.files(path='~/PregFCS',pattern='.fcs',full.names=TRUE)
source('~/Pregnancy/Tools/Process_PP.R')
################################################

MN=readRDS('Processed/MN_Pregnancy.rds')
ToUse=readRDS('Processed/ToUse_Pregnancy.rds')
Meta_Preg=readRDS('Processed/Meta_Pregnancy.rds')

#################
#clustering
################
startT=Sys.time()
Build=runRepMetaclust(50,50,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)
endT=Sys.time()
CT=endT-startT
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

startT=Sys.time()
AUCs=runClassif(FuncDF,Meta_Preg$Class,40,IterNumClus,0.5,1,as.character(Meta_Preg$Subject),35)
endT=Sys.time()
fullT=endT-startT
print('done with classification')
stop('')
#############
#Visualization
##############
CellMat_Pregy=readRDS('Processed/CellMat_Preg.rds')
print('starting visualization')
saveDir='OutDir/Pregnancy_Viz'

Layout_Pregnancy=Rtsne(CelMat_Preg)$Y
Atlas=vizAtlas(CellMat_Preg,Build,Meta_Preg$Class,ToUse,SampsToUse=NULL,35,Layout_Pregnancy,saveDir)

