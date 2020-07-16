vizFunctionMaps=function(Layout,FuncFeat,MN,funcMarker,Resp,outDir){
	#purpose: color clusters generated through VoPo clustering by difference in functional marker expression between groups
	#Inputs:
		#Layout: The layout of clusters in 2D returned by vizClusterPhenotype.R
		#FuncFeat: The matrix of function-based features created by getFunctionalFeature.R	
		#MN: the vector of comprehensible marker names corresponding to the channels of the FCS files
		#funcMarker: A vector of indices corresponding to the functional markers.
			#Ex: For example if CREB and NFKB were in channels 2 and 3 and we wanted to look at differences for each then funcMarker=c(2,3)
		#Resp: The response variables corresponding to the rownames of FuncFeat (and also FileName input to VoPo). Assumes this is a binary classification problem with 2 classes
		#outDir: The directory where we should save this plots

	#Output: This function will write the function-based plots to your outDir. There should be length(funcMarker) plots, corresponding to each marker of interest

library('ggplot2')
library('reshape2')
library('viridis')


#get possible response values
vals=unique(Resp)
Inds0=which(Resp==vals[1])

MarkerNames=MN[funcMarker]
MarkNames=colnames(FuncFeat)

for(f in 1:length(funcMarker)){

fInds=grep(MarkerNames[f],MarkNames)
fInds=sort(fInds)

SubF=FuncFeat[,fInds]

wRes=c()
for(c2 in 1:ncol(SubF)){
	AVals=SubF[Inds0,c2]
	BVals=SubF[-Inds0,c2]

wRes=c(wRes,wilcox.test(AVals,BVals,alternative=c('two.sided'))$p.value)
} #for c2

wRes=log(wRes,10)
tooLow=which(wRes<(-6))
wRes[tooLow]=-6
wRes[is.na(wRes)]=0

MarkName2=MarkerNames[f]

NetDF=data.frame(Layout[,1],Layout[,2],wRes)
 names(NetDF)=c('PC1','PC2','log10pVal')
  p=ggplot(NetDF, aes(PC1,PC2,color=log10pVal))+
     geom_point() +
      scale_color_viridis(option='B',direction=-1,limits=c(-6,0))+ theme_bw()+ theme(text = element_text(size=20))+ggtitle(MarkName2)
  p=p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

FName=paste(outDir,'/',MarkName2,'.jpg',sep='')
 ggsave(FName,p,width=7,height=7)
} #for f





}
