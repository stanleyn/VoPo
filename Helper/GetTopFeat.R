GetTopFeat=function(X,NumKeep,NN){

library('igraph')

#Step 1 is to create kNN graph
NNGraph=kNN(X,NN)

#Step 2 is to build similarity matrix
S=matrix(0,nrow=nrow(NNGraph),ncol=nrow(NNGraph))
for(ss in 1:nrow(S)){
	Edges=which(NNGraph[ss,]==1)
	for(j in 1:length(Edges)){
		if(ss<Edges[j]){
		x1=X[ss,]
		x2=X[Edges[j],]
		diff=norm(as.matrix(x1-x2),type='F')
		diff=(diff^2)/5
		S[ss,Edges[j]]=exp(-diff)
		S[Edges[j],ss]=exp(-diff)
		}
	}
}

#Step 3: Compute Laplacian
GetLaplace=Laplacian(S)
L=GetLaplace[[2]]
D=GetLaplace[[1]]

#Generate matrix of ones
Ones=as.matrix(rep(1,nrow(S)),ncol=1)

#Step 4: Compute score for each feature
FScore=c()
for(r in 1:ncol(X)){
	#get their feature vector
	fr=as.matrix(X[,r],ncol=1)
	Num=t(fr)%*%D%*%Ones
	Denom=t(Ones)%*%D%*%Ones
	Subtract=c(Num/Denom)*Ones
	TfR=fr-Subtract

	#Compute Laplacian Score
	Num2=t(TfR)%*%L%*%TfR
	Denom2=t(TfR)%*%D%*%TfR
	Lr=Num2/Denom2
	FScore=c(FScore,Lr)

}

#get indices for the top scoring features
TopFeat=order(FScore,decreasing=TRUE)[1:NumKeep]
#print(TopFeat)
Out=X[,TopFeat]
Out
} ##function end 

#########################################
#Helper Functions
###########################################
#Purpose: Helper functions for feature selection
#Date: July 19

###########
#kNN graph#
###########
kNN=function(Mat,NN){
k=NN
dist_mat <- as.matrix(dist(Mat, method = "euclidean", upper = TRUE, diag=TRUE))
nrst <- lapply(1:nrow(dist_mat), function(i) k.nearest.neighbors(i, dist_mat, k = k))
w <- matrix(nrow = dim(dist_mat), ncol=dim(dist_mat)) ## all NA right now
w[is.na(w)] <- 0 ## populate with 0
for(i in 1:length(nrst)) for(j in nrst[[i]]) w[i,j] = 1
Adj2=w
Net=graph.adjacency(Adj2,mode='undirected')
FinalAdj=get.adjacency(Net,type='both',sparse=FALSE)
diag(FinalAdj)=0
FinalAdj	 
}

##########
#Laplacian#
##########
Laplacian=function(Adj){
Ds=rowSums(Adj)
D=matrix(0,nrow=nrow(Adj),ncol=nrow(Adj))
diag(D)=Ds
L=D-Adj
Out=list()
Out[[1]]=D
Out[[2]]=L
Out
}
