#Purpose: Regression with random forest with multiview feature selection and robust CV pipeline
#Clean version for public use. :)
#Date: April 13, 2020
#Dependency: please also install package MLmetrics
#Author: Natalie Stanley stanleyn@stanford.edu
runRegressionTask=function(FuncDF,Y,FPV,IterNumClus,propTrain,numPerm,sampID,numCore){

#dependencies
library('randomForest')
library('foreach')
library('doParallel')
library('matrixStats')
library('ROCR')
source('Helper/MultiviewFS.R')
source('Helper/GetTopFeat.R')
library('FastKNN')
library('MLmetrics')
#inputs:
	#FuncDF: The sample x feature per iteration matrix returned by B2
	#Y: the continuous response vector
	#FPV: The number of features to use per view. 
	#IterNumClus: The number of clusters per iteration
	#propTrain: the proportion of the data to train with
	#sampID: Your input will either be 0 or a vector of subject IDs for each sample. Input a vector of subjectIDs for each sample if every subject has more than 1 FCS file.
	#numCore: The number of cores to use for parallelization
#Output:
	#StoreMSEs: the vector of MSEs corresponding to each cross validation trial

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

StoreMSEs=c()

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
trainGet=sample(UniqueSamp,length(UniqueSamp)*propTrain)

#get all sample IDs that is an element of this
relSampInds=which(is.element(sampID,trainGet))
train=relSampInds

test=seq(length(Y))[-train]
ansVec=rep(NA,length(Y))
mod=randomForest(X2[train,],Y[train])
predVec=predict(mod,X2[test,])
ansVec[test]=predVec
return(list(ansVec=ansVec))
} #RF iter loop

stopCluster(cl)

ansMat=RFIter$ansVec
ans=rowMedians(ansMat,na.rm=TRUE)

#calculate the mean squared error between true and actual predictions
r=MSE(ans,Y)
StoreMSEs=c(StoreMSEs,r)
} #for ss
StoreMSEs
} #function end
