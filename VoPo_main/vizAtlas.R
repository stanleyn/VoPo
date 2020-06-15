vizAtlas=function(CellMat,Build,Y,ToUse,SampsToUse=NULL,numCore,layout,outDir){
	library('foreach')
	library('doParallel')
	library('miscTools')
	library('ggplot2')
	library('reshape2')
	library('viridis')
	library('pROC')
#description:
#This function is meant to directly correspond to the repeated metaclustering analysis
#inputs:
#CellMat: The cell x marker matrix. A limited subsample of cells across files
#Build is the data structure returned from repeated metaclustering (for script runRepMetaclustering.R
#Y: is the response variable corresponnding to the order of FileNames given to runRepMetaclustering
#ToUse: the indicies of the functional and phenotypic markers you what to use in your 2d layout
#UseSamps: all-> use all samples the calculate statistics specify-> use specific indices to calculate statistics
#SampsToUse: The indices of samples to use if UseSamps='specify'       		
#numCore: the number of cores to use in the parallelization
#layout: This should be the 2D coordinates corresponding to the cells in CellMat 
#outDir: the path to the directory where the output plots should go

	#Get the cell medians
	ClusMeds=Build$ClusMed[,ToUse]

	#Get the frequencies
	FreqVals=Build$FProp

	#get IterNumClus
	IterNumClus=Build$IterNumClus

	if(length(SampsToUse)>1){
	FreqSelect=FreqVals[SampsToUse,]
	SelectResponse=Y[SampsToUse]
	}
	else{
	FreqSelect=FreqVals
	SelectResponse=Y
	}
	
	#prepare for pvalue calculation
	U1=unique(Y)[1]
	aInds=which(SelectResponse==U1)
	pVal=c()

	S=length(IterNumClus)

	#comopute p value for each cluster
	for(i in 1:sum(IterNumClus[1:S])){
		aVals=FreqSelect[aInds,i]
		bVals=FreqSelect[-aInds,i]
		ww=wilcox.test(aVals,bVals)$p.value
		pVal=c(pVal,ww)
	}

	#start building the visualization over all bootstraps
	cl=makeCluster(numCore)
	registerDoParallel(cl)
	PointMat=foreach(s=1:S,.combine='rbind',.packages=c('FNN')) %dopar% {

	OutMat=rep(0,nrow(CellMat))

	#grab the relevant rows for metacluster centers
	if(s==1){
		start=1
	}
	else{
	start=sum(IterNumClus[1:(s-1)])+1 
	}          
	end=sum(IterNumClus[1:s])
	
	subMed=ClusMeds[start:end,]

	#now calculate the distance for each center to 
	numCl=IterNumClus[s]
	
	for(d in 1:IterNumClus[s]){
		
		#FNN implementation
		distValsTemp=knnx.dist(CellMat,t(subMed[d,]),k=nrow(CellMat))
		forOrder=knnx.index(CellMat,t(subMed[d,]),k=nrow(CellMat))
		distVals=rep(0,length(distValsTemp))
		distVals[forOrder]=distValsTemp
		distVals=exp(-5*distVals)
		OutMat=rbind(OutMat,distVals)
	} #end d

	OutMat=OutMat[-1,]
	
	#now do calculation for linear combination of pvalues
	subPVal=pVal[start:end]
	pointVal=c()

	for(i in 1:ncol(OutMat)){
	vec=OutMat[,i]
	tempVal=sum(vec*subPVal)/sum(vec)
	pointVal=c(pointVal,tempVal)
		}
return(pointVal)
}
stopCluster(cl)
print('starting plot')
#############
##Plotting###
############
setwd(outDir)

#first make plots of phenotypic markers
clFeat=CellMat
for(i in 1:ncol(clFeat)){
NetDF=data.frame(layout[,1],layout[,2],clFeat[,i])
names(NetDF)=c('PC1','PC2','Marker')
FName=paste(colnames(clFeat)[i],'.jpg',sep='')
p=ggplot(NetDF, aes(PC1,PC2,color=Marker))+geom_point()+ scale_color_viridis()+ theme_bw()+ theme(text = element_text(size=20))+xlab('dim1')+ylab('dim 2')+ggtitle(colnames(clFeat)[i])+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(FName,p,width=7,height=7)
   }

#color points by their computed significance score log-scale
pointVal=colMeans(log(PointMat,10))

 NetDF=data.frame(layout[,1],layout[,2],pointVal)
 names(NetDF)=c('PC1','PC2','pointVal')
  p=ggplot(NetDF, aes(PC1,PC2,color=pointVal))+geom_point(size=.1)+ scale_color_viridis(option='B',direction=-1)+ theme_bw()+ theme(text = element_text(size=20))
   p=p+theme(axis.title.x = element_blank(),
	       axis.title.y = element_blank())
   p=p+theme(axis.line = element_line(colour = "black"),
	         panel.grid.major = element_blank(),
		     panel.grid.minor = element_blank(),
		     panel.border = element_blank(),
		         panel.background = element_blank())+theme(legend.title=element_blank())
     p=p+theme(legend.position='bottom')+xlab('dim1')+ylab('dim2')
   ggsave('pval.jpg',p,width=7,height=7)


PointMat
}
