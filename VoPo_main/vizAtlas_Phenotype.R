vizAtlas_Phenotype=function(CellMat,tRes,outDir){
library('viridis')
library('ggplot2')
library('reshape2')
setwd(outDir)
#first make plots of phenotypic markers
clFeat=CellMat
layout=tRes
for(i in 1:ncol(clFeat)){
NetDF=data.frame(layout[,1],layout[,2],clFeat[,i])
names(NetDF)=c('PC1','PC2','Marker')
FName=paste(colnames(clFeat)[i],'.jpg',sep='')
p=ggplot(NetDF, aes(PC1,PC2,color=Marker))+geom_point()+ scale_color_viridis()+ theme_bw()+ theme(text = element_text(size=20))+xlab('dim1')+ylab('dim 2')+ggtitle(colnames(clFeat)[i])+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(FName,p,width=7,height=7)
   }
}
