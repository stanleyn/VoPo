
# Overview and Instructions

## Updates and Contact Information

* Date Updated: June 15, 2020
* Prepared By: Natalie Stanley (stanleyn@stanford.edu)
* Code Tested in: R version 3.4.4 

## Purpose

We introduce VoPo which enables end-to-end bioinformatics analysis of single-cell data. Here we provide code for general usage and to reproduce the results in our paper `VoPo Leverages Cellular Heterogeneity for Predictive Modeling of Single-Cell Data`. The (4) tasks outlined in this readme enable the following:

1) **Task 1:Clustering and Feature Engineering**: Run VoPo clustering on some FCS files

2) **Task 2: Classification and Reproduce Fig 2D**: Feature Use VoPo features for classification tasks in each of the 3 datasets. Also reproduce the results in Figure 2D. Generating distributions of classification accuracies from single vs. repeated metaclustering solutions (Figure 2D.). 

3) **Task 3: Differentiation Score Visualization (Figs 2a-c)**: Generate comprehensive single-cell visualizations for each of the 3 datasets (Figure 2A-C). Uses processed data from running 50 iterations of the repeated metaclustering algorithm. 

Alternatively, in (4) we have provided an example for re-running all clustering results from scratch. However this will require downloading FCS files from flow repository and modifying the paths in the scripts to the data, accordingly.  

4) **Task 4: Re-run repeated metaclustering from scratch on the 3 datasets**: Re-running the repeated metaclustering strategy from scratch along with classification and visualization. Note that this task requires downloading FCS files from each dataset. 

5) **Task 5: Visualize differentiation score on single cells based on frequency or function-based features**: Examples on how to create single-cell visualizations with points colored by differentiation score (Fig 2 a-c) for frequency features and for function features. 

## Dependencies

Dependencies: We use a number of R packages. Please make sure that you have these packages installed. `flowCore`, `FNN`, `foreach`, `doParallel`, `iterators`, `plyr`, `randomForest`, `matrixStats`,`ROCR`,`FastKNN`, `miscTools`, `ggplot2`, `reshape2`, `viridis`, `pROC`, `igraph`.

## Installation Instructions

* We assume that you have R installed. :) We used R version 3.4.4 to produce the results in this paper. R and the above packages are the only things that you need to install and we expect this to take minuts on a standard desktop computer.

* You can clone this git repository by 

```Bash
git clone https://github.com/stanleyn/VoPo
```
* You cloned this git repository into some place, `*YourPath*`. Once you are in R, please change your working directory in R so that you are in this folder, as all paths are relative to here.

```R
setwd('*YourPath*')
```
## License
This software is licensed under Apache License, Version 2.0 (https://www.apache.org/licenses/LICENSE-2.0). 

# Demos and Reproduction of Results

## Task 1: Run VoPo Clustering and Extract Features

* Here is a script for a quick demo for how to run VoPo clustering on FCS files. Once you have run this script, use the output to perform classification and visualization tasks (see Tasks 2-4). 

* This demo should take less than 5 minutes to run on a standard computer. We used 5 cores in this example.

* Please see `Demo_VoPo` script for a description of the inputs, etc.

```R
source('Demo_Data/Demo_VoPo.R')
```

* This script created a VoPo object called `Build` that can be used for further tasks.

* If you want to extract frequency-related features from VoPo clustering object, `Build`, and the vector of file names you gave to VoPo `FNames` in the above Demo_VoPo script, you can do this with the following:

```R
source('VoPo_main/getFrequencyFeature.R')
FrequencyFeatures=getFrequencyFeature(Build,FNames)
```

This is the data matrix you can use for classification tasks (like we did in the paper). However, in the below example I show you how you can also get function-related features.

* You can also extract function based features. You will input `Build` (the VoPo clustering object), `FNames`, which are the filenames in the order you gave them to VoPo (see Demo_VoPo script) and `FInds`, which are the indices of functional markers in your FCS files

```R
source('VoPo_main/getFunctionalFeature')
#Here are the indices of the functional markers
FInds=readRDS('Processed/FI_HSR.rds')
FunctionalFeatures=getFunctionalFeature(Build,FNames,FInds)
```

## Task 2: Generate Distributions of Classification Accuracies for Single vs. Repeated Metaclustering Solutions (Fig 2D.)

* For each dataset, we will show you how to generate a distribution of classification accuracies from the VoPo engineered features. For each dataset, the first example referring to the script in the `Examples` directory shows how you can build a model based on the extracted features.
* We will also generate the boxplots (baseline to VoPo distribution comparison in figure 2D and figures with appear in the 'OutDir' directory.

### Hip Surgery Recovery Dataset (HSR)

Get a distribution of VoPo classification accuracies for HSR dataset
```R
source('Examples/Class_Surgery_Example.R')
```

The resulting vector of classification accuracies is `ClAcc`

We can also create the boxplots from figure 2D in the HSR dataset.

```R
source('PaperFigures/Classification/Class_Surgery.R')
```
You can now find your boxplots in OutDir as SurgDist.pdf

### Normal Term Pregnancy Dataset (NTP)

Get a distribution of VoPo classification accuracies for the NTP dataset

```R
source('Examples/Class_Surgery_Example.R')
```

The resulting vector of classification accuracies is `ClAcc`.

We can also create the boxplots for Figure 2D in the NTP dataset.

```R
source('PaperFigures/Classification/Class_Pregnancy.R')
```
You can now find your boxplots in OutDir as Preg_Dist.pdf

### Longitudinal Stroke Recovery Dataset (LSR)

Get a distribution of VoPo classification accuracies in the LSR dataset.

```R
source('Examples/Class_Stroke_Example.R')
```

The resulting vector of classification accuracies is `ClAcc`

We can also generate boxplots for the HSR dataset from Figure 2D

```R
source('PaperFigures/Classification/Class_Stroke.R')
```
You can now find your boxplots in OutDir as Stroke_Dist.pdf

## Task 3: Comprehensive Single-Cell Visualizations for each Dataset (Fig 2A-C.)

* Results for each dataset will be within their respective folder in OutDir. 
* There is a plot for each surface marker showing its expression across cells (ex. CD3.jpg shows CD3 expression)
* pval.jpg plot colors the cells by their differentiation scores computed by the algorithm. 

### Hip Surgery Recovery Dataset (HSR)

```R
source('PaperFigures/Visualization/SurgeryViz.R')
```
You can now find plots for all markers and differentiation scores in OutDir/Surgery_Viz

### Normal Term Pregnancy Dataset (NTP)

```R
source('PaperFigures/Visualization/PregnancyViz.R')
```
You can now find plots for all markers and differentiation scores in OutDir/Pregnancy_Viz

### Longitudinal Stroke Recovery Dataset (LSR)

```R
source('PaperFigures/Visualization/StrokeViz.R')
```

You can now find plots for all markers and differentiation scores in OutDir/Stroke_Viz

## Task 4: Re-run clustering from scratch. Use that clustering result to generate classification results and visualizations (Fig 2D.)

* You can run this part of the code if you have access to the raw FCS files. You will need to go in and change the path to the FCS files depending on where you download them. The place where you need to change the path is clearly marked :)

* These scripts will re-run the classification and visualization pipelines

### Hip Surgery Recovery Dataset (HSR)

Make sure that the order of your Files matches these names

```R
Meta_Surgery=readRDS('Processed/Meta_Surgery')
Meta_Surgery$FileNames
```

**Change the path in the below script to wherever your FCS files are!!!!!**

```R
source('PaperFigures/RunClustering/SurgeryRe.R')
```

* You can find the visualization results in OutDir/Surgery_Viz

* You can find the AUCs from 30 runs of the cross validation pipeline stores in the vector `AUCs`

### Normal Term Pregnancy Dataset (NTP)

Make sure that the order of your Files matches these names

```R
Meta_Preg=readRDS('Processed/Meta_Pregnancy')
Meta_Preg$FileNames
```

**Change the path in the below script to wherever your FCS files are!!!!**

```R
source('PaperFigures/RunClustering/PregnancyRe.R')
```

* You can find the visualization results in OutDir/Pregnancy_Viz

* You can find the AUCs from 30 runs of the cross validation pipeline stores in the vector `AUCs`

### Longitudinal Stroke Recovery Dataset (LSR)

Make sure that the order of your Files matches these names

```R
Meta_Stroke=readRDS('Processed/Meta_Stroke')
Meta_Stroke$FileNames
```

**Change the path in the below script to wherever your FCS files are!!!!**

```R
source('PaperFigures/RunClustering/StrokeRe.R')
```

You can find the visualization results in OutDir/Stroke_Viz

You can find the AUCs from 30 runs of the cross validation pipeline stores in the vector `AUCs`

## Task 5: Create single-cell visualization with each cell computed by differentiation score. 

* Here we will show examples of how to color by differentation score for Frequency Based Features and Functional Features.
* I will tell you the functions to use here. **Examples coming soon!** Inputs are similar to what you have seen above.
* A differentiation score plot based on frequency will ultimately create a single plot.
* Using functional marker features, we will get one differentiation score plot for each functional marker.

Instructions for making differentiation score plots for both frequency and function

* Step 1: Select a large subset of cells across all of your FCS files. Use function, 

```R
source('VoPo_main/SampleCells.R') 
```

* Step 2: Make a frequency map. Each cell will be colored by differentiation score according to frequency. Use the function,

```R
source('VoPo_main/vizFrequencyMap.R')
```

* Step 3: Make functional maps. For each functional marker in your panel, each cell will be colored according to differentiation score based on the expression of that particular functional marker. Note that as one of the functional arguments you will specify the directory to write these plots to. They will all go to one place. Use the function,

```R
source('VoPo_main/vizFunctionMaps.R')
```

