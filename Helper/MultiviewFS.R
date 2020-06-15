#Purpose: The purpose of this code is to do multiview embedding for prediction
#In this code we assume we are using the median data for pregnancy
#Input:
	#NumIt: the number of iterations to use
	#DataMat: the data matrix we will use to concat
	#ClusNum: the number of clusters
	#FPV: The number of features per view
	#NN: the number of nearest neighbors to use in the graph

#step 1: make the concatenated PCA matrix
MultiviewFS=function(NumIt,DataMat,ClusNum,FPV,NN){
ConcatPCA=rep(0,nrow(DataMat))

for(i in 1:NumIt){
if(i==1){
low=1
}

else{low=sum(ClusNum[1:(i-1)])+1
}

high=sum(ClusNum[1:i])

IterMat=DataMat[,low:high]

#prevent from asking for more columns than we have
FPV2=min(ncol(IterMat),FPV)

Y1=GetTopFeat(IterMat,FPV2,NN)

if(ClusNum[i]==0){
ConcatPCA=ConcatPCA
}

else{
ConcatPCA=cbind(ConcatPCA,Y1)}

} #for

ConcatPCA=ConcatPCA[,-1]

##############################
ConcatPCA
}
