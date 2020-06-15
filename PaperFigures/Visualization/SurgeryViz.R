#Date: October 11, 2019
#Purpose: Create the phenotypic and significance visualization for surgery dataset

source('VoPo_main/vizAtlas.R')
source('VoPo_main/vizAtlas_FS.R')

library('FastKNN')

Build_Surgery=readRDS('Processed/Build_Surgery_NoGrans.rds')
CellMat_Surgery=readRDS('Processed/CellMat_Surgery.rds')
Layout_Surgery=readRDS('Processed/layouts/kehlet_Y.rds')
ToUse_Surgery=readRDS('Processed/ToUse_Surgery.rds')
Meta_Surgery=readRDS('Processed/Meta_Surgery.rds')
SampInds_Surg=readRDS('Processed/relSampInds_Surg.rds')

saveDir='OutDir/Surgery_Viz'
Atlas=vizAtlas(CellMat_Surgery,Build_Surgery,Meta_Surgery$Class,ToUse_Surgery,SampsToUse=SampInds_Surg,35,Layout_Surgery,saveDir)
stop('')

######
#Please contact authors for other kinds of analyses
######

#using largeVis coordinates
#lay=readRDS('~/Clean_BClust/LargeVis_April7/LVEmbed_Kehlet')
#lay=t(lay)

#lay=cbind(CellMat_Surgery[,8],CellMat_Surgery[,3])
#Atlas=vizAtlas(CellMat_Surgery,Build_Surgery,Meta_Surgery$Class,ToUse_Surgery,SampsToUse=SampInds_Surg,35,lay,saveDir)
#stop('')

#for fold change
#Fold=CellViz_FoldChange(CellMat_Surgery,Build_Surgery,Meta_Surgery$Class,ToUse_Surgery,SampInds_Surg,50)


#Atlas=vizAtlas_FS(CellMat_Surgery,Build_Surgery,Meta_Surgery$Class,ToUse_Surgery,SampsToUse=SampInds_Surg,35,Layout_Surgery,saveDir,40)


