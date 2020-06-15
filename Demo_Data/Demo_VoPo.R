#Purpose: This is the steps needed to run bootsrapped BL2 on Kehlet
#Date Created: July 25, 2018
######################################################

#dependencies
library('flowCore')
library('foreach')
library('doParallel')
library('iterators')
library('Rtsne')

source('VoPo_main/runRepMetaclust.R')

########################################
##Step 1: Put in filenames, info etc. #
########################################

#Put the path to your FCS files
FileNames=list.files(path='Demo_Data',pattern='.fcs',full.names=TRUE)

#These are the human read-able antibody names corresponding to the columns of each FCS file
MN=readRDS('Processed/MN_Surgery.rds')

#Note tha we have provided the pre-processed version of this here, but you can extract them as follows
#frame=read.FCS(FileNames[1]
#MN=pData(parameters(frame))[,2]

#To use gives the indices, corresponding to the columns of FCS files that you want to use for clustering
ToUse=readRDS('Processed/ToUse_Surgery.rds')

Meta_Surgery=readRDS('Processed/Meta_Surgery.rds')

##############################################
#run clustering
###############################################

#inputs in the order that they appear
	#50: The number of metaclustering iterations
	#50: The number of metaclusters per iteration
	#FileNames: The vector of filenames with the full path defined above
	#doCPF='specify' means that we are going to specify the # of clusters per file
	#numCPF: The number of clusters per FCS file
	#MN: The vector of marker names that correspond to known antibody names that we either loaded or defind above
	#ToUse: The vector of indices corresponding to the columns we want to use for clustering
	#5: We are running this with 5 cores.

Build=runRepMetaclust(50,50,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,5)
print('clustering done!')

