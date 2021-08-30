vizAtlas_Freq_Directional=function(CellMat,Build,FreqMat,Y,ToUse,SampsToUse=NULL,numCore,layout,outDir,numIter){
	library('foreach')
	library('doParallel')
	library('miscTools')
	library('ggplot2')
	library('reshape2')
	library('viridis')
	library('pROC')
	library('scales')

#description:
#This function is meant to directly correspond to the repeated metaclustering analysis
#does a directional analysis
#inputs:
#CellMat: The cell x marker matrix. A limited subsample of cells across files
#Build is the data structure returned from repeated metaclustering (for script runRepMetaclustering.R
#FreqMat: the matrix of frequency features you made using getFrequencyFeature.R
#Y: is the response variable corresponnding to the order of FileNames given to runRepMetaclustering
#ToUse: the indicies of the functional and phenotypic markers you what to use in your 2d layout
#SampsToUse: The indices of samples to use. If you want to use all, keep default NULL       	
#numCore: the number of cores to use in the parallelization
#layout: This should be the 2D coordinates corresponding to the cells in CellMat 
#outDirFunc: the path to the directory where the output plots should go
#outDitPhen: the path to the directory where phenotype plots should go
#FuncNames: Names of functional markers you would like to plot
#numIter: The number of clustering iterations to use

	#Get the cell medians
	ClusMeds=Build$ClusMed[,ToUse]

	IterNumClus=Build$IterNumClus

	
	if(length(SampsToUse)>1){
	FreqSelect=FreqMat[SampsToUse,]
	SelectResponse=Y[SampsToUse]
	}
	else{
	FreqSelect=FreqMat
	SelectResponse=Y
	}
	
	#prepare for pvalue calculation
	U1=min(Y)[1]
	print(U1)
	
	aInds=which(SelectResponse==U1)

	S=length(IterNumClus)

	wRes=c()
        for(c2 in 1:ncol(FreqSelect)){
	
	AVals=FreqSelect[aInds,c2]
	BVals=FreqSelect[-aInds,c2]
	
	testVal=wilcox.test(AVals,BVals,alternative=c('two.sided'))$p.value
	
	MeanA=mean(AVals,na.rm=TRUE) #mean for A class
	MeanB=mean(BVals,na.rm=TRUE) #mean for B class

	#take log of p val and replace NaN with 0
	testVal=log(testVal,10)
	testVal[is.nan(testVal)]=0
	testVal[testVal<(-5)]=-4

	#subtract mean to figure out sign
	SubVal=MeanB-MeanA

	#correct sign by direction
	getSign=sign(SubVal)

	#multiply by -1 so it will be positive values
	testVal=-1*testVal

	#correct by sign
	testVal=testVal*getSign

	wRes=c(wRes,testVal)
	} #for c2

	pVal=wRes
	nanInds=which(is.nan(pVal))
	pVal[nanInds]=0
	print(min(pVal))
	print(max(pVal))
	#start building the visualization over all bootstraps
	cl=makeCluster(numCore)
	registerDoParallel(cl)
	
	#specify s so that it only goes over the number of iterations you actually want
	S=numIter

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
FName=paste('frequency','.jpg',sep='')
#color points by their computed significance score log-scale
pointVal=colMeans(PointMat,10)
print(max(pointVal))
print(min(pointVal))

 NetDF=data.frame(layout[,1],layout[,2],pointVal)
 names(NetDF)=c('PC1','PC2','pointVal')
  p=ggplot(NetDF, aes(PC1,PC2,color=pointVal))+geom_point(size=.1)+scale_color_gradient2(low = muted("blue"), mid = "gray95",high = muted("red"), midpoint = 0, space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "colour",limits=c(-4,4))+ theme_bw()+ theme(text = element_text(size=14))
   p=p+theme(axis.title.x = element_blank(),
	       axis.title.y = element_blank())
   p=p+theme(axis.line = element_line(colour = "black"),
	         panel.grid.major = element_blank(),
		     panel.grid.minor = element_blank(),
		     panel.border = element_blank(),
		         panel.background = element_blank())+theme(legend.title=element_blank())
     p=p+theme(legend.position='bottom')+xlab('dim1')+ylab('dim2')
   ggsave(FName,p,width=7,height=7)


}
