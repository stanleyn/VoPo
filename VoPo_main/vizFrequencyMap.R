vizFrequencyMap=function(FreqFeat,Layout,Resp,outDir){
	#purpose: visualize frequency differences in each of the VoPo identified clusters
	#Inputs:
		#FreqFeat: Frequency-based features obtained using getFrequencyFeature.R
		#Layout: The Layout returned from vizClusterPhenotype.R
		#Resp: The response variables that correspond to the rows of FreqFeat. Assumes binary response variable
		#outDir: The path where this plot should be written to
	#Outputs:
		#A plot of frequency differences in VoPo identified clusters will be saved to outDir


library('ggplot2')
library('reshape2')
library('viridis')


#get possible response values
vals=unique(Resp)
Inds0=which(Resp==vals[1])

PVVec=c()
 for(f in 1:ncol(FreqFeat)){
 	stat=wilcox.test(FreqFeat[Inds0,f],FreqFeat[-Inds0,f])$p.value
 	PVVec=c(PVVec,stat)
 }

#clean up 0s
zeroInds=which(PVVec==0)
PVVec[zeroInds]=0.0000001
PVVec=log(PVVec,10)

 NetDF=data.frame(Layout[,1],Layout[,2],PVVec)
 names(NetDF)=c('PC1','PC2','log10PVal')
 p=ggplot(NetDF, aes(PC1,PC2,color=log10PVal))+
  geom_point() +
  scale_color_viridis(option='B',direction=-1)+ theme_bw()+ theme(text = element_text(size=20))+ggtitle('Frequency Diffs')
p=p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
FName=paste(outDir,'/','Freq.pdf',sep='')
ggsave(FName,p,width=7,height=7)


}
