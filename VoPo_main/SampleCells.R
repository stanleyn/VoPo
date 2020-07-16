SampleCells=function(FileNames,MN,ToUse,NumCells,NumCellsFinal){
library('doParallel')
library('foreach')
library('flowCore')
#Date: Updated June 25, 2020
#Description: This is the function for building our iter
#Inputs (only 1 so far):
	#FileNames: The full path name for the file names
	#MN: The names of the markers (measured CyTOF parameters) in order they appea in files
	#ToUse: The indices of functional markers to use for clustering
	#NumCells: The number of cells to take from each sample
	#NumCellsFinal: The ultimate number of cells (e.g. 30,000)
#Outputs
	#$PatMark: Functional marker exp in each iteration cluster per patient
	#$PatProp: Frequency data--> prop of cells in each cluster for each sample
	
#set up the parallelization
cl=makeCluster(20)
registerDoParallel(cl)

#creating a custom combine for S
custom_s=function(L1,L2){
PatMark=cbind(L1$PatMark,L2$PatMark)
PatProp=cbind(L1$PatProp,L2$PatProp)
ConsFeatByMark=rbind(L1$ConsFeatByMark,L2$ConsFeatByMark)
return(list(PatMark=PatMark,PatProp=PatProp,ConsFeatByMark=ConsFeatByMark))
#return(list(PatMark=PatMark,PatProp=PatProp))
}

custom_i=function(L1,L2){
	SubDM=rbind(L1$SubDM,L2$SubDM)
	AllMarker=rbind(L1$AllMarker,L2$AllMarker)
	return(list(SubDM=SubDM,AllMarker=AllMarker))
}

##############################################
#Step 1: Define data matrix used in clustering
##############################################
DM_Build=foreach(i=FileNames,.combine=custom_i) %do% {

#process file
cytoftrans=arcsinhTransform(transformationId='cytofTransform',a=0,b=(1/5),c=0)
frame=read.FCS(i)

#transformation and matrix building
translist=transformList(colnames(exprs(frame)),cytoftrans)
frame=transform(frame,translist)
X=exprs(frame)

#sample NumCells cells from the file
NCU=min(NumCells,nrow(X))
SampInds=sample(1:nrow(X),NCU,replace=FALSE)
AllMarker=X[SampInds,]
SubDM=X[SampInds,ToUse]
colnames(SubDM)=MN[ToUse]
colnames(AllMarker)=MN
return(list(SubDM=SubDM,AllMarker=AllMarker))
} #foreach i

AllMarker=DM_Build$AllMarker

#return only the sub data matrix
ClDM=DM_Build$SubDM

#now subsample NumCellsFinal number of cells
SampInds=sample(1:nrow(ClDM),NumCellsFinal,replace=FALSE)
print(dim(ClDM))
ClDM=ClDM[SampInds,]
ClDM
}
