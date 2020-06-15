vizAtlas_FS=function(CellMat,Build,Y,ToUse,SampsToUse=NULL,numCore,layout,outDir,FPV){
	library('foreach')
	library('doParallel')
	library('miscTools')
	library('ggplot2')
	library('reshape2')
	library('viridis')
	library('pROC')
	library('Rtsne')
source('~/Clean_BClust/General/MultiviewFS.R')
source('~/Clean_BClust/General/GetTopFeat.R')

#description:
#Important note: this function is to visualize which features are included
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

	NewName=paste(1:ncol(FreqVals),colnames(FreqVals),sep='--')
	colnames(FreqVals)=NewName

	##stuff that we need for PCA version
	#Layout=prcomp(scale(ClusMeds))$x[,1:2]
	#create cor network
	#Cor=t(FreqVals)%*%FreqVals
	#Layout=Rtsne(Cor,check_duplicates=FALSE)$Y

	print('starting')

#	cl=makeCluster(numCore)
#	registerDoParallel(cl)
	
#	PointMat=foreach(s=1:50,.combine='rbind',.packages=c('FNN')) %dopar% {
#	library('FastKNN')
#	source('~/Clean_BClust/General/MultiviewFS.R')
#	source('~/Clean_BClust/General/GetTopFeat.R')
	
	#sample half of the data
	#numSamp=floor(0.5*nrow(FreqVals))
	#sampVals=sample(1:nrow(FreqVals),numSamp,replace=FALSE)
	#X1=MultiviewFS(50,FreqVals[sampVals,],Build$IterNumClus,FPV,10)
 		
	#relInds=which(is.element(colnames(FreqVals),colnames(X1)))
	#pVal=rep(0,ncol(FreqVals))
	#pVal[relInds]=1
	#return(pVal)
	#}

#print(PointMat)
#stopCluster(cl)

#PointMat=colSums(PointMat)/nrow(PointMat)

#do the plotting
# NetDF=data.frame(Layout[,1],Layout[,2],PointMat)
#print(dim(NetDF)) 
#names(NetDF)=c('PC1','PC2','pointVal')
#  p=ggplot(NetDF, aes(PC1,PC2,color=pointVal))+geom_point(size=1)+ scale_color_gradient(low='#E1FA72', high='#F46FEE')+ theme_bw()+ theme(text = element_text(size=20))
#   p=p+theme(axis.title.x = element_blank(),
#       axis.title.y = element_blank())
#   p=p+theme(axis.line = element_line(colour = "black"),
#	         panel.grid.major = element_blank(),
#		     panel.grid.minor = element_blank(),
#		     panel.border = element_blank(),
#		         panel.background = element_blank())+theme(legend.title=element_blank())
#     p=p+theme(legend.position='bottom')+xlab('dim1')+ylab('dim2')
#   ggsave('FS_include.jpg',p,width=7,height=7)

#g=ggplot_build(p)
#coloring=g$data[[1]]['colour']
#  colors=coloring[,1]

#library('igraph')
#library('FastKNN')

#make KNN
#k=20
#    dist_mat <- as.matrix(dist(ClusMeds, method = "euclidean", upper = TRUE, diag=TRUE))
#   nrst <- lapply(1:nrow(dist_mat), function(i) k.nearest.neighbors(i, dist_mat, k = k))
#       w <- matrix(nrow = dim(dist_mat), ncol=dim(dist_mat)) ## all NA right now
#       w[is.na(w)] <- 0 ## populate with 0
#           for(i in 1:length(nrst)) for(j in nrst[[i]]) w[i,j] = 1

       # #  #create the network object
 #         Adj2=w
#	     Net=graph.adjacency(Adj2,mode='undirected')	  
#	Net=simplify(Net)
#	  layout=layout_with_fr(Net,coord=Layout) #layout with fr is the best so far
#	   #quartz()
#	   postscript('~/FS_include.eps',width=6,height=6)
#	   ColVec=colors
#	   plot(Net,vertex.label=NA,vertex.size=3,edge.color='gray88',vertex.color=ColVec,layout=Layout,vertex.frame.color=NA)
#dev.off()

############### #end potential stochastic use of the code############


	S=length(Build$IterNumClus)
	IterNumClus=Build$IterNumClus

	#feature selection part
	X1=MultiviewFS(50,FreqVals,Build$IterNumClus,FPV,10)

	#figure out which original column names were retained
	relInds=which(is.element(colnames(FreqVals),colnames(X1)))
	pVal=rep(0,ncol(FreqVals))
	pVal[relInds]=1
#	print(sum(pVal))

#	print(relInds)

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

#for(i in 1:ncol(clFeat)){
#NetDF=data.frame(layout[,1],layout[,2],clFeat[,i])
#names(NetDF)=c('PC1','PC2','Marker')
#FName=paste(colnames(clFeat)[i],'.jpg',sep='')
#p=ggplot(NetDF, aes(PC1,PC2,color=Marker))+geom_point()+ scale_color_viridis()+ theme_bw()+ theme(text = element_text(size=20))+xlab('dim1')+ylab('dim 2')+ggtitle(colnames(clFeat)[i])+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
#ggsave(FName,p,width=7,height=7)
 #  }

#color points by their computed significance score log-scale
pointVal=colMeans(PointMat)

NetDF=data.frame(layout[,1],layout[,2],pointVal)
 names(NetDF)=c('PC1','PC2','pointVal')
  p=ggplot(NetDF, aes(PC1,PC2,color=pointVal))+geom_point(size=.1)+ scale_color_gradient(low='#E1FA72', high='#F46FEE')+ theme_bw()+ theme(text = element_text(size=20))
   p=p+theme(axis.title.x = element_blank(),
	       axis.title.y = element_blank())
   p=p+theme(axis.line = element_line(colour = "black"),
	         panel.grid.major = element_blank(),
		     panel.grid.minor = element_blank(),
		     panel.border = element_blank(),
		         panel.background = element_blank())+theme(legend.title=element_blank())
     p=p+theme(legend.position='bottom')+xlab('dim1')+ylab('dim2')
   ggsave('FS_include.jpg',p,width=7,height=7)


#PointMat
}
#PointMat
#}
