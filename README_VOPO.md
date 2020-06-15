
# Overview

## Updates and Contact Information

* Date Updated: June 15, 2020
* Prepared By: Natalie Stanley (stanleyn@stanford.edu)
* Code Tested in: R version 3.4.4 (2018-03-15)

## Purpose

We introduce VoPo which enables end-to-end bioinformatics analysis of single-cell mass cytometry data. Here we provide code to reproduce the following results. 

1) Running VoPo clustering on a few FCS files

2) Generating distributions of classification accuracies from single vs. repeated metaclustering solutions (Figure 2D.). Uses processed data from running 50 iterations of the repeated clustering algorithm.

3) Generating comprehensive single-cell visualizations for each of the 3 clinical datasets (Figure 2A-C). Uses processed data from running 50 iterations of the repeated clustering algorithm. 

Alternatively, in (4) we have provided an example for re-running all clustering results from scratch. However this will require downloading FCS files from flow repository and modifying the paths in the scripts to the data, accordingly.  

4) Re-running the repeated metaclustering strategy from scratch along with classification and visualization. Note that this task requires downloading FCS files from each dataset. 

## Dependencies

Dependencies: We use a number of R packages. Please make sure that you have these packages installed. `flowCore`, `foreach`, `doParallel`, `iterators`, `plyr`, `randomForest`, `matrixStats`,`ROCR`,`FastKNN`, `miscTools`, `ggplot2`, `reshape2`, `viridis`, `pROC`, `igraph`.

# Installation Instructions

* We assume that you have R installed. :) We used R version 3.4.4 to produce the results in this paper. R and the above packages are the only things that you need to install and we expect this to take minuts on a standard desktop computer.

* You can clone this git repository by 

```Bash
git clone https://github.com/stanleyn/VoPo
```
* You cloned this git repository into some place, `*YourPath*`. Once you are in R, please change your working directory in R so that you are in this folder

```R
setwd('*YourPath')
```
Demos:

*************************
#########################
#Task 1: Run VoPo clustering on some FCS files
*************************
########################

Here is a script for a quick demo for how to run VoPo clustering on FCS files. Once you have run this script, use the output to perform classification and visualization tasks (see Tasks 2-4). 

This demo should take less than 5 minutes to run on a standard computer. We used 5 cores in this example.

>setwd('YourPath/Reproduce')
>source('Demo_Data/Demo_VoPo.R')
>MyOutput=Build

#This creates a list called Build, which is the processed data that you can use for further tasks. 

**************************
##########################
#Task 2: Generate Distributions of Classification Accuracies for Single vs. Repeated Metaclustering Solutions (Fig 2D.)
#########################
***************************
Assuming you are in the Reproduce directory, we will show how to create the boxplots shown in Fig2D for each dataset

##################################
#Hip Surgery Recovery Dataset (HSR)
##################################
>setwd('YourPath/Reproduce')
>source('PaperFigures/Classification/Class_Surgery.R')

You can now find your boxplots in OutDir as SurgDist.pdf

##################################
#Normal Term Pregnancy Dataset (NTP)
###################################
>setwd('YourPath/Reproduce')
>source('PaperFigures/Classification/Class_Pregnancy.R')

You can now find your boxplots in OutDir as Preg_Dist.pdf

###################################
#Longitudinal Stroke Recovery Dataset (LSR)
###################################
>setwd('YourPath/Reproduce')
>source('PaperFigures/Classification/Class_Stroke.R')

You can now find your boxplots in OutDir as Stroke_Dist.pdf

**************************
##########################
#Task 3: Comprehensive Single-Cell Visualizations for each Dataset(Fig 2D.)
##########################
***************************

-Results for each dataset will be within their respective folder in OutDir. 
-There is a plot for each surface marker showing its expression across cells (ex. CD3.jpg shows CD3 expression)
-pval.jpg plot colors the cells by their differentiation scores computed by the algorithm. 

##################################
#Hip Surgery Recovery Dataset (HSR)
##################################
>setwd('YourPath/Reproduce')
>source('PaperFigures/Visualization/SurgeryViz.R')

You can now find plots for all markers and differentiation scores in OutDir/Surgery_Viz

##################################
#Normal Term Pregnancy Dataset (NTP)
###################################
>setwd('YourPath/Reproduce')
>source('PaperFigures/Visualization/PregnancyViz.R')

You can now find plots for all markers and differentiation scores in OutDir/Pregnancy_Viz

###################################
#Longitudinal Stroke Recovery Dataset (LSR)
###################################
>setwd('YourPath/Reproduce')
>source('PaperFigures/Visualization/StrokeViz.R')

You can now find plots for all markers and differentiation scores in OutDir/Stroke_Viz

**************************
##########################
#Task 4: Re-run clustering from scratch. Use that clustering result to generate classification results and visualizations (Fig 2D.)
#########################
***************************

-You can run this part of the code if you have access to the raw FCS files. You will need to go in and change the path to the FCS files depending on where you download them. The place where you need to change the path is clearly marked :)
-These scripts will re-run the classification and visualization pipelines

##################################
#Hip Surgery Recovery Dataset (HSR)
##################################

>setwd('YourPath/Reproduce')

-Make sure that the order of your Files matches these names
>Meta_Surgery=readRDS('Processed/Meta_Surgery')
>Meta_Surgery$FileNames

-Change the path in the below script to wherever your FCS files are!!!!!

>source('PaperFigures/RunClustering/SurgeryRe.R')

-You can find the visualization results in OutDir/Surgery_Viz
-You can find the AUCs from 30 runs of the cross validation pipeline stores in the vector AUCs

>AUCs 

##################################
#Normal Term Pregnancy Dataset (NTP)
###################################
>setwd('YourPath/Reproduce')

-Make sure that the order of your Files matches these names
>Meta_Preg=readRDS('Processed/Meta_Pregnancy')
>Meta_Preg$FileNames

-Change the path in the below script to wherever your FCS files are!!!!

>source('PaperFigures/RunClustering/PregnancyRe.R')

-You can find the visualization results in OutDir/Pregnancy_Viz
-You can find the AUCs from 30 runs of the cross validation pipeline stores in the vector AUCs

>AUCs 

###################################
#Longitudinal Stroke Recovery Dataset (LSR)
###################################
>setwd('YourPath/Reproduce')

-Make sure that the order of your Files matches these names
>Meta_Stroke=readRDS('Processed/Meta_Stroke')
>Meta_Stroke$FileNames

-Change the path in the below script to wherever your FCS files are!!!!

>source('PaperFigures/RunClustering/StrokeRe.R')

-You can find the visualization results in OutDir/Stroke_Viz
-You can find the AUCs from 30 runs of the cross validation pipeline stores in the vector AUCs

>AUCs

#############################################
#############################################
#License
#############################################
#############################################

This software is licensed under Apache License, Version 2.0 (https://www.apache.org/licenses/LICENSE-2.0). 
