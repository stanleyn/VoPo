getFrequencyFeature=function(Build,FileNames){

#purpose: turn VoPo output into a nice matrix of frequency based features
	#If there are K clusters per VoPo iteration and I VoPo iterations, there will be I x K total features 

#inputs:
	#Build: the output returned from VoPo clustering (from running runRepMetaClust.R)
	#FileNames: The vector of fileNames that you input to runRepMetaClust.R. This is for rownames of your matrix
#output:
	#FreqMat: the matrix of frequency based features from VoPo. Dimensions will be length(FileNames) x total features



FreqDF=Build$FProp

#modify the column names to make sure they are all unique
ForC=1:ncol(FreqDF)
NewName=paste(colnames(FreqDF),ForC,sep='_')
colnames(FreqDF)=NewName
rownames(FreqDF)=FileNames

FreqDF
}
