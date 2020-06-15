runRepMetaclust=function(S,K,FileNames,doCPF=c('auto','specify'),numCPF=NULL,MN,ToUse,numCore){
#Date: Updated October 7, 2019
#Description: Metaclustering code to partition cells into k metaclusters
#Inputs:
	#S: The number of repeated metaclustering iterations to do (e.g. 50)
	#K: The number of metaclusters (e.g. 50)
	#FileNames: The vector of FCS file names (e.g. c('File1.fcs','File2.fcs')). Use the full path to the files
	#doCPF: whether or not to choose the number of clusters per file automatically (doCPF='auto') or specify (doCPF='specify')
	#numCPF: number of clusters per file. Leave as NULL if doCPF='auto'
	#MN: A vector of marker names (measured CyTOF markers) corresponding to columns of FCS files (e.g. c('CD4','CD8','FoxP3'))
	#ToUse: the indices corresponding to the markers to be used for clustering (e.g. c(1,2,3,4))

library('flowCore')
library('foreach')
library('doParallel')
library('iterators')
library('plyr')

cl=makeCluster(numCore)
registerDoParallel(cl)

#custom combine for within-file clustering
custom=function(LL1,LL2){
Assn=c(LL1$Assn,LL2$Assn)
Centers=rbind(LL1$Centers,LL2$Centers)
UClus=c(LL1$UClus,LL2$UClus)
NumCells=c(LL1$NumCells,LL2$NumCells)
XSample=rbind(LL1$XSample,LL2$XSample)
return(list(XSample=XSample,Assn=Assn,Centers=Centers,UClus=UClus,NumCells=NumCells))
}

DM_Cluster=foreach(i=FileNames,.combine=custom,.packages=c('flowCore','plyr')) %dopar% {

###################
####process file###
###################
cytoftrans=arcsinhTransform(transformationId='cytofTransform',a=0,b=(1/5),c=0)
frame=read.FCS(i)

####for applying transformation
translist=transformList(colnames(exprs(frame)),cytoftrans)
frame=transform(frame,translist)
X=exprs(frame)
X=X[,1:length(MN)]
colnames(X)=MN

############################
######within file clustering#
#############################

X2=scale(X) 
X2=X2[,ToUse] 

if(doCPF=='auto'){
NumCl=floor(sqrt((nrow(X)/2)))
}

else{NumCl=min(numCPF,nrow(X2)-1)}

ClRes=kmeans(X2,centers=NumCl)

#collect results
Assn=ClRes$cluster
Centers=ClRes$centers
UClus=nrow(Centers)
NumCells=nrow(X)
Assn2=as.factor(Assn)
X=data.frame(X,Assn2)
names(X)[ncol(X)]='Assn2'
XSample=ddply(X,.(Assn2),numcolwise(median))
XSample=as.matrix(XSample[,-1])
return(list(XSample=XSample,Assn=Assn,Centers=Centers,UClus=UClus,NumCells=NumCells))
}

stopCluster(cl)

################################
#begin metaclustering, feature construction, etc.
###########################

cl=makeCluster(numCore)
registerDoParallel(cl)

custom_s=function(L1,L2){
FMed=cbind(L1$FMed,L2$FMed)
FProp=cbind(L1$FProp,L2$FProp)
IterNumClus=c(L1$IterNumClus,L2$IterNumClus)
ClusMed=rbind(L1$ClusMed,L2$ClusMed)
return(list(FMed=FMed,FProp=FProp,IterNumClus=IterNumClus,ClusMed=ClusMed))
}

SecondDM=DM_Cluster$Centers
SecondDM2=scale(SecondDM)
XSampFull=DM_Cluster$XSample

print('starting metaclustering')

BootMeta=foreach(s=1:S,.combine='custom_s',.packages=c('foreach','flowCore','plyr')) %dopar%{

NumClus2=K

#actual metaclustering
MetaClust=kmeans(SecondDM2,centers=NumClus2)
MetaAssn=MetaClust$cluster
MetaCenter=MetaClust$centers

#calculate the median marker expression for each cluster
toAssn=as.factor(MetaAssn)
XX=data.frame(XSampFull,toAssn)
names(XX)[ncol(XX)]='Assn'
ClusMed=ddply(XX,.(toAssn),numcolwise(median))
ClusMed=as.matrix(ClusMed[,-1])
IterNumClus=NumClus2          

########################################################
#mapping metacluster labels back to single cellls
#####################################################

#Indiviudal file cluster centers were assigned to metaclusters
#Individual file clusters are comprised of cells in each file]
#We will map the metacluster labels back to cells
UClus=DM_Cluster$UClus
Assn=DM_Cluster$Assn
NumCells=DM_Cluster$NumCells

#Assn keeps track of within-file cluster assignment for each cell across files
#We will update sample-specific entries of Assn to metacluster labels
#Index bookkeeping
old=0
cellOld=0
startEndMat=matrix(0,nrow=length(UClus),ncol=2)
for(i in 1:length(UClus)){
Start=old+1
End=Start+(UClus[i]-1)
old=End
startCell=cellOld+1
endCell=startCell+(NumCells[i]-1)
cellOld=endCell
startEndMat[i,1]=startCell
startEndMat[i,2]=endCell
RelAssn=Assn[startCell:endCell]
clIndexes=Start:End
Convert=clIndexes[RelAssn]
Assn[startCell:endCell]=Convert
}

#CellMetaLab is the converted version ot cell-to-metacluster labels
CellMetaLab=rep(0,length(Assn))

c_i=function(L1,L2){
calcMed=rbind(L1$calcMed,L2$calcMed)
propVec=rbind(L1$propVec,L2$propVec)
CellMetaLab=cbind(L1$CellMetaLab,L2$CellMetaLab)
return(list(calcMed=calcMed,propVec=propVec,CellMetaLab=CellMetaLab))
}

c_j=function(L1,L2){
calcMed=c(L1$calcMed,L2$calcMed)
propVec=c(L1$propVec,L2$propVec)
return(list(calcMed=calcMed,propVec=propVec))
}

Base2PFeat=foreach(i=1:length(FileNames),.combine=c_i,.packages='foreach') %do% {
startCell=startEndMat[i,1]
endCell=startEndMat[i,2]
GetAssn=Assn[startCell:endCell]
WorkWith=MetaAssn[GetAssn]
CellMetaLab[startCell:endCell]=WorkWith
ClsMeta=c(1:max(MetaAssn))

###########################################################
#begin calculating frequency and functional marker features
###########################################################

foreach(j=ClsMeta,.combine=c_j) %do% {
relInds=which(WorkWith==j)
stage2Assn=GetAssn[relInds]
tabStage2Assn=table(stage2Assn)

#keeps track of how many cells in each file were assigned to each individual cluster
Weights=t(as.matrix(tabStage2Assn))
IndexWeight=as.numeric(names(tabStage2Assn))

if(length(unique(stage2Assn))==0){
calcMed=rep(0,ncol(XSampFull))
}

else if(length(unique(stage2Assn))==1){
calcMed=as.matrix(XSampFull[IndexWeight,])
calcMed=calcMed[,1]
}

else{
#get weighted mean
WeightMed=Weights%*%as.matrix(XSampFull[IndexWeight,])
calcMed=WeightMed/sum(Weights)
calcMed=calcMed[1,]
}

#get names for medians
Name4Med=paste(j,colnames(XSampFull),sep='_')
names(calcMed)=Name4Med

#calculate frequencies
propVec=length(relInds)/NumCells[i] 
names(propVec)=paste(j,'--prop',sep='')
return(list(calcMed=calcMed,propVec=propVec))
	} #j

} #i
FMed=Base2PFeat$calcMed
FProp=Base2PFeat$propVec
return(list(FMed=FMed,FProp=FProp,IterNumClus=IterNumClus,ClusMed=ClusMed))
} #s

stopCluster(cl)
BootMeta
}

