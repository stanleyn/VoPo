#Purpose: This is the steps needed to run bootsrapped BL2 on Kehlet
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
##Step 1: Put in filenames, info etc. #
########################################

#Modify here!!! Put the path to your FCS files
FileNames=list.files(path='~/FCS',pattern='.fcs',full.names=TRUE)
###################################

MN=readRDS('Processed/MN_Surgery.rds')
ToUse=readRDS('Processed/ToUse_Surgery.rds')
Meta_Surgery=readRDS('Processed/Meta_Surgery.rds')

##############################################
#run the entire pipeline
###############################################

############
#clustering
############
ClStart=Sys.time()
Build=runRepMetaclust(50,50,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)
ClEnd=Sys.time()
ClTime=ClEnd-ClStart
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

Start=Sys.time()
AUCs=runClassif(FuncDF,Meta_Surgery$Class,40,IterNumClus,0.5,1,as.character(Meta_Surgery$Subject),35)
endT=Sys.time()
Class=endT-Start
print('classification done!')

stop('')
##############
#Visualization
##############
SampInds_Surg=readRDS('Processed/relSampInds_Surg.rds')
CellMat_Surgery=readRDS('Processed/CellMat_Surgery.rds')

#run tSNE
Layout_Surgery=Rtsne(Cell_Mat_Surgery)

saveDir='OutDir/Surgery_Viz'
Atlas=vizAtlas(CellMat_Surgery,Build,Meta_Surgery$Class,ToUse,SampsToUse=SampInds_Surg,35,Layout_Surgery,saveDir)
