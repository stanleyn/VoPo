#Purpose: Classification with multiview feature selection and robust CV pipeline
#Clean version for public use. :)
#Date: July 4, 2019
#Author: Natalie Stanley stanleyn@stanford.edu
runClassif=function(FuncDF,Y,FPV,IterNumClus,propTrain,numPerm,sampID,numCore){

#dependencies
library('randomForest')
library('foreach')
library('doParallel')
library('matrixStats')
library('ROCR')
source('Helper/MultiviewFS.R')
source('Helper/GetTopFeat.R')
library('FastKNN')
#inputs:
	#FuncDF: The sample x feature per iteration matrix returned by B2
	#Y: the binary response vector
	#FPV: The number of features to use per view. 
	#IterNumClus: The number of clusters per iteration
	#propTrain: the proportion of the data to train with
	#sampID: Your input will either be 0 or a vector of subject IDs for each sample. Input a vector of subjectIDs for each sample if every subject has more than 1 FCS file.
	#numCore: The number of cores to use for parallelization
#Output:
	#StoreAUCs: the vector of AUCs

#######################################################

####################################
######Parameter Specification#######
####################################

#Parameters for feature selection
N=10

#Use the maximal number of clustering iterations provided
m=length(IterNumClus)

#specify the sampID#
if(length(sampID)==1){
sampID=c(1:nrow(FuncDF))
}

else{sampID=sampID}
####################################
####################################

StoreAUCs=c()

for(ss in 1:numPerm){

####################################################
#####Random Forest construction#####
###################################################

#set up parallelization
cl=makeCluster(numCore)
registerDoParallel(cl)

rfCustom=function(L1,L2){
ansVec=cbind(L1$ansVec,L2$ansVec)
modImp=rbind(L1$modImp,L2$modImp)
return(list(ansVec=ansVec,modImp=modImp))
}

#number of CV iterations to do
itSeq=1:500

X2=MultiviewFS(m,FuncDF,IterNumClus,FPV,N)

RFIter=foreach(s=itSeq,.combine=rfCustom,.packages='randomForest') %dopar% {

UniqueSamp=unique(sampID)

repeat{
trainGet=sample(UniqueSamp,length(UniqueSamp)*propTrain)
print(trainGet)
#get all sample IDs that is an element of this
relSampInds=which(is.element(sampID,trainGet))

train=relSampInds
if(length(unique(Y[train]))==2)
break
}

test=seq(length(Y))[-train]
ansVec=rep(NA,length(Y))
mod=randomForest(X2[train,],Y[train])
predVec=predict(mod,X2[test,],type='prob')[,2]
ansVec[test]=predVec
return(list(ansVec=ansVec))
} #RF iter loop

stopCluster(cl)

ansMat=RFIter$ansVec
ans=rowMedians(ansMat,na.rm=TRUE)

#calculate AUC
pred=prediction(ans,Y)
r1=performance(pred,measure='auc')
r=r1@y.values[[1]]
print('trial #')
print(ss)
print(r)
StoreAUCs=c(StoreAUCs,r)
} #for ss

StoreAUCs
} #function end
