FS_returnFeat=function(FuncDF,IterNumClus,FPV){

#purpose: return the indices of features selected by Laplacian scoring feature selection
#inputs:
	#FuncDF: matrix of frequency-based features from repeated metaclustering
	#IterNumClus: Number of clusters per iteration from VoPo result
	#FPV: the number of clusters to select per solution
#outputs:
	#vector of indices for the features that should be used

source('Helper/MultiviewFS.R')
source('Helper/GetTopFeat.R')
library('FastKNN')

colnames(FuncDF)=paste(colnames(FuncDF),c(1:ncol(FuncDF)),sep='_')

X2=MultiviewFS(length(IterNumClus),FuncDF,IterNumClus,FPV,10)
relInds=which(is.element(colnames(FuncDF),colnames(X2)))

relInds
}

