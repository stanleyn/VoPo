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

#Modify here!!! Put the path to your FCS files
FileNames=list.files(path='Demo_Data',pattern='.fcs',full.names=TRUE)
###################################

MN=readRDS('Processed/MN_Surgery.rds')
ToUse=readRDS('Processed/ToUse_Surgery.rds')
Meta_Surgery=readRDS('Processed/Meta_Surgery.rds')

##############################################
#run the entire pipeline
###############################################

############
#clustering
############
ClStart=Sys.time()
Build=runRepMetaclust(50,50,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,5)
ClEnd=Sys.time()
ClTime=ClEnd-ClStart
print('clustering done!')

