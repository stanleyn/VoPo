library('flowCore')
library('foreach')
library('doParallel')
library('iterators')
library('Rtsne')

source('VoPo_main/runRepMetaclust.R')
source('VoPo_main/runClassif.R')
source('VoPo_main/vizAtlas.R')

###############################
##Step 1: Put in filename, info, etc.
##################################

#Modify here!!! Put the path to your FCS files
#FileNames=list.files('~/Stroke_FCS')
source('~/Stroke/Process_Stroke_June.R')
FileNames=FNames

###############################################

MN=readRDS('Processed/MN_Stroke.rds')
ToUse=readRDS('Processed/ToUse_Stroke.rds')
Meta_Stroke=readRDS('Processed/Meta_Stroke.rds')

###############################
#Clustering
###############################
startT=Sys.time()
Build=runRepMetaclust(50,50,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)
endT=Sys.time()
CT=endT-startT
print('done with clustering')

####################
#Classification
################
IterNumClus=Build$IterNumClus
FuncDF=Build[[2]]
colnames(FuncDF)=1:ncol(FuncDF)

startT=Sys.time()
AUCs=runClassif(FuncDF,factor(Meta_Stroke$Class),40,IterNumClus,0.5,1,as.character(Meta_Stroke$Subject),35)
endT=Sys.time()
ClassT=endT-startT

print('done with classification')
stop('')
#############################
#Visualization
#############################
CellMat_Stroke=readRDS('Processed/CellMat_Stroke.rds')
saveDir='OutDir/Stroke_Viz'
Layout_Stroke=Rtsne(CellMat_Stroke)$Y
Atlas=vizAtlas(CellMat_Stroke,Build,factor(Meta_Stroke$Class),ToUse,SampsToUse=NULL,35,Layout_Stroke,saveDir)
