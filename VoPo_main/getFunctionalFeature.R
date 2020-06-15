getFunctionalFeature=function(Build,FileNames,FuncInds){

#purpose: turns VoPo output into a nice matrix of function based features
	#If you have F functional markers, K clusters per iteration and, I VoPo iterations, there are F*K*I total features

#inputs:
	#Build: The output from VoPo clustering (runRepMetaClust.R)
	#FileNames: The vector of fileNames that you input to runRepMetaClust.R. This is for rownames of your matrix
	#FuncInds: A vector with the indices of the functional markers
		#forexample: if columns 1, 5, 7 were functional markers FuncInds=c(1,5,7)


#output: a matrix of size length(FileNames) x total number of features

FuncDF=Build$FMed
CN=as.matrix(colnames(Build$ClusMed)[FuncInds])

ARes=apply(CN,1,function(x) which(grepl(x,colnames(FuncDF))))
UInds=sort(unique(unlist(ARes)))
FuncDF=FuncDF[,UInds]

#modify the column names to make sure they are all unique
ForC=1:ncol(FuncDF)
NewName=paste(colnames(FuncDF),ForC,sep='_')
colnames(FuncDF)=NewName
rownames(FuncDF)=FileNames

FuncDF
}
