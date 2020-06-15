#Date: October 11, 2019
#Purpose: Create the phenotypic and significance visualization for pregnancy dataset

source('VoPo_main/vizAtlas.R')
source('VoPo_main/vizAtlas_FS.R')
library('FastKNN')

#To Update!!!::::read in info for pregnancy dataset
Build_Preg=readRDS('Processed/Build_Preg_NoGrans.rds')

CellMat_Preg=readRDS('Processed/CellMat_Pregnancy.rds')
Layout_Preg=readRDS('Processed/layouts/pregnancy_Y.rds')
ToUse_Preg=readRDS('Processed/ToUse_Pregnancy.rds')
Meta_Preg=readRDS('Processed/Meta_Pregnancy.rds')
		  
saveDir='OutDir/Pregnancy_Viz'

Atlas=vizAtlas(CellMat_Preg,Build_Preg,Meta_Preg$Class,ToUse_Preg,SampsToUse=NULL,35,Layout_Preg,saveDir)
stop('')

##################################################
#for additional analyses: please contact authors
################################################

#fold change
#Fold=CellViz_FoldChange(CellMat_Preg,Build_Preg,Meta_Preg$Class,ToUse_Preg,1:62,50)
#stop('')

#using largeVis coordinates
#lay=readRDS('~/Clean_BClust/LargeVis_April7/LVE0708_Preg')
#lay=t(lay$coord)
#lay=cbind(CellMat_Preg[,6],CellMat_Preg[,25])

#Atlas=vizAtlas(CellMat_Preg,Build_Preg,Meta_Preg$Class,ToUse_Preg,SampsToUse=NULL,35,lay,saveDir)
#stop('')

#Atlas=vizAtlas_FS(CellMat_Preg,Build_Preg,Meta_Preg$Class,ToUse_Preg,SampsToUse=NULL,35,Layout_Preg,saveDir,40)

#stop('')

