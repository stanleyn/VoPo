B2_adaptive=function(S,K,FileNames,MN,ToUse,Cl_Out){
#Date: Updated June 15, 2020
#Description: This is our adaptive metaclustering
#Inputs:
	#S: The number of iterations
	#K, the number of clusters in the metaclustering step -> 0 is auto
	#FileNames: The full path for the gile names
	#MN: The names of markers (measured CyTOF parameters)
	#ToUse: the indices of phenotypic markers used for clustering
	#Cl_Out: put in the cluster centers from already obtained clustering results
library('flowCore')
library('foreach')
library('doParallel')
library('iterators')
library('plyr')
library('FNN')

cl=makeCluster(35)
registerDoParallel(cl)

#####################################################
#########Data specific stuff:surgery#####
####################################################

#creating a custom combine for S
custom_s=function(L1,L2){
FMed=cbind(L1$FMed,L2$FMed)
FProp=cbind(L1$FProp,L2$FProp)
return(list(FMed=FMed,FProp=FProp))
}

#creating a custom combine for building DM of cluster centers
custom=function(LL1,LL2){
Assn=c(LL1$Assn,LL2$Assn)
Centers=rbind(LL1$Centers,LL2$Centers)
UClus=c(LL1$UClus,LL2$UClus)
NumCells=c(LL1$NumCells,LL2$NumCells)
XSample=rbind(LL1$XSample,LL2$XSample)
return(list(XSample=XSample,Assn=Assn,Centers=Centers,UClus=UClus,NumCells=NumCells))
}
DM_Cluster=foreach(i=FileNames,.combine=custom,.packages=c('flowCore','plyr')) %dopar% {
print(i)

####process file
cytoftrans=arcsinhTransform(transformationId='cytofTransform',a=0,b=(1/5),c=0)
frame=read.FCS(i)

####for applying transformation
translist=transformList(colnames(exprs(frame)),cytoftrans)
frame=transform(frame,translist)
X=exprs(frame)
##keep this only for the quirk of the pregnancy data
X=X[,1:length(MN)]
###################
colnames(X)=MN
print(dim(X))
############################
######within file clustering#
#############################

#figure out how many clusters ####change to autostart
#NumCl=floor(sqrt((nrow(X)/2)))

#number of clusters should be generalizeable even if there is low number
NumCl=min(1000,nrow(X)-1)

#NumCl=2000
#perform clustering
X2=scale(X)  #this is the scaled version
X2=X2[,ToUse] 
ClRes=kmeans(X2,centers=NumCl)

######putt out the data structures we need
Assn=ClRes$cluster
Centers=ClRes$centers

######track how many cells and centers we have per file
UClus=nrow(Centers)
NumCells=nrow(X)

##create a data frame of the data and 
Assn2=as.factor(Assn)
X=data.frame(X,Assn2)
names(X)[ncol(X)]='Assn2'
##create the per cluster XSample
XSample=ddply(X,.(Assn2),numcolwise(median))
XSample=as.matrix(XSample[,-1])

#Determine the start and end inds
return(list(XSample=XSample,Assn=Assn,Centers=Centers,UClus=UClus,NumCells=NumCells))
}
#stop the cluster from the original file step
stopCluster(cl)

#start a new one for the repeated clustering step
cl=makeCluster(35)
registerDoParallel(cl)

custom_s=function(L1,L2){
FMed=cbind(L1$FMed,L2$FMed)
FProp=cbind(L1$FProp,L2$FProp)
IterNumClus=c(L1$IterNumClus,L2$IterNumClus)
ClusMed=rbind(L1$ClusMed,L2$ClusMed)
return(list(FMed=FMed,FProp=FProp,IterNumClus=IterNumClus,ClusMed=ClusMed))
}

#get stuff we need for the metaclustering
SecondDM=DM_Cluster$Centers
SecondDM2=scale(SecondDM)

#if(K==0){
#NumClus2=floor(sqrt((nrow(SecondDM)/2)))}

#else if(K=='samp'){
#NumClus2=sample(5:30,1)
#}

#else{NumClus2=K}

###also grab the X sample so that we have the full centers
XSampFull=DM_Cluster$XSample
print('starting metaclustering')
BootMeta=foreach(s=1:S,.combine='custom_s',.packages=c('FNN','foreach','flowCore','plyr')) %dopar% {

if(K==0){
NumClus2=floor(sqrt((nrow(SecondDM)/2)))}


else if(K=='samp'){
NumClus2=sample(5:30,1)
}

else{NumClus2=K}

#this is the over clustering step
MetaClust=kmeans(SecondDM2,centers=NumClus2)

#get thing we are going to match it to
startInd=(K*(s-1))+1
endInd=K*s

MetaCenter=Cl_Out[startInd:endInd,ToUse]

#match Second DM2 to MetaCenter
#forOrder=knnx.index(MetaCenter,SecondDM2,k=1)
#MetaAssn=c(forOrder)

clRes=kmeans(SecondDM2,centers=NumClus2)
MetaAssn=clRes$cluster

#make sure there is at least 1 thing assigned to each k
repreVec=length(unique(MetaAssn))
repeat{
for(kk in 1:K){
Is=which(MetaAssn==kk)
if(length(Is)==0){
sampInds=sample(1:nrow(SecondDM2),1)
MetaAssn[sampInds]=kk
}
}
repreVec=length(unique(MetaAssn))
if(repreVec==K){
break
}
}

#calculate the cluster x median matrix
toAssn=as.factor(MetaAssn)
XX=data.frame(XSampFull,toAssn)
names(XX)[ncol(XX)]='Assn'
ClusMed=ddply(XX,.(toAssn),numcolwise(median))
ClusMed=as.matrix(ClusMed[,-1])


IterNumClus=NumClus2          #getNumClus

#map back to the original cells
UClus=DM_Cluster$UClus
Assn=DM_Cluster$Assn
NumCells=DM_Cluster$NumCells

##initialize the previous index 
old=0
cellOld=0

#create a start end mat
startEndMat=matrix(0,nrow=length(UClus),ncol=2)

for(i in 1:length(UClus)){
#figure out indices for cluster number 
Start=old+1
End=Start+(UClus[i]-1)
old=End

#get the corresponding cell indices to that file
startCell=cellOld+1
endCell=startCell+(NumCells[i]-1)
cellOld=endCell

#update the start end Mat
startEndMat[i,1]=startCell
startEndMat[i,2]=endCell

#pull out the relevant Assn indices
RelAssn=Assn[startCell:endCell]

#create vector of indices of clusters
clIndexes=Start:End

#convert to these new values
Convert=clIndexes[RelAssn]

#replace with new values
Assn[startCell:endCell]=Convert
}

CellMetaLab=rep(0,length(Assn))

#done mapping: start building feature matrices
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
#get the correct indices for the start-- end cells
startCell=startEndMat[i,1]
endCell=startEndMat[i,2]
#get the relevant Assningments for this
GetAssn=Assn[startCell:endCell]
#convert these to the metaclusters
WorkWith=MetaAssn[GetAssn]
#this is the cell to metadata label
CellMetaLab[startCell:endCell]=WorkWith

#enumerate cluster #s from metaassn
ClsMeta=c(1:max(MetaAssn))

foreach(j=ClsMeta,.combine=c_j) %dopar% {
#get the indices of the pulled indexes in class j
relInds=which(WorkWith==j)
#map back to their original assignments in stage 2
stage2Assn=GetAssn[relInds]
#tablulate the results
tabStage2Assn=table(stage2Assn)
#get the weights for each class
Weights=t(as.matrix(tabStage2Assn))
#get the indices for the above weights
IndexWeight=as.numeric(names(tabStage2Assn))

if(length(unique(stage2Assn))==0){
calcMed=rep(0,ncol(XSampFull))
}

else if(length(unique(stage2Assn))==1){
calcMed=as.matrix(XSampFull[IndexWeight,])
calcMed=calcMed[,1]
}

else{
#get weighted median
WeightMed=Weights%*%as.matrix(XSampFull[IndexWeight,])
calcMed=WeightMed/sum(Weights)
calcMed=calcMed[1,]
}

##get names for median
Name4Med=paste(j,colnames(XSampFull),sep='_')
names(calcMed)=Name4Med

####proportion calculation
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

