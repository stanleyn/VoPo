#Date: October 11, 2019
#Purpose: Create the phenotypic and significance visualization for stroke dataset

source('VoPo_main/vizAtlas.R')
source('VoPo_main/vizAtlas_FS.R')
library('FastKNN')

#read in info for stroke dataset
Build_Stroke=readRDS('Processed/Build_Stroke_NoGrans.rds')

CellMat_Stroke=readRDS('Processed/CellMat_Stroke.rds')
Layout_Stroke=readRDS('Processed/layouts/stroke_Y.rds')
ToUse_Stroke=readRDS('Processed/ToUse_Stroke.rds')
Meta_Stroke=readRDS('Processed/Meta_Stroke.rds')

saveDir='OutDir/Stroke_Viz'
Atlas=vizAtlas(CellMat_Stroke,Build_Stroke,Meta_Stroke$Class,ToUse_Stroke,SampsToUse=NULL,35,Layout_Stroke,saveDir)
stop('')

########
#for additional analyses: please contact authors
#######

#Fold=CellViz_FoldChange(CellMat_Stroke,Build_Stroke,Meta_Stroke$Class,ToUse_Stroke,1:40,50)

#lay=readRDS('~/Clean_BClust/LargeVis_April7/LVEmbed_Stroke')
#lay=t(lay$coord)

#lay=cbind(CellMat_Stroke[,8],CellMat_Stroke[,6])
#Atlas=vizAtlas(CellMat_Stroke,Build_Stroke,Meta_Stroke$Class,ToUse_Stroke,SampsToUse=NULL,35,lay,saveDir)

#stop('')
#Atlas=vizAtlas_FS(CellMat_Stroke,Build_Stroke,Meta_Stroke$Class,ToUse_Stroke,SampsToUse=NULL,35,Layout_Stroke,saveDir,40)

