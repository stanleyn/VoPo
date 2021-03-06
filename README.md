
# Overview and Instructions

## Updates and Contact Information

* Date Updated: June 15, 2020
* Prepared By: Natalie Stanley (stanleyn@stanford.edu)
* Code Tested in: R version 3.4.4 
* **Important Note: make sure you have git lfs installed! Some of the data to run the tutorials are stored in git lfs, and without this, the files will be corrupted** (https://git-lfs.github.com/)

## Purpose

We introduce VoPo which enables end-to-end bioinformatics analysis of single-cell data. Here we provide code for general usage and to reproduce the results in our paper `VoPo Leverages Cellular Heterogeneity for Predictive Modeling of Single-Cell Data`. **You can focus on tasks 1-4 for instructons on general usage: clustering, feature extraction, and visualization (e.g. coloring cells by differentiation score)**. Tasks 1-4 outlined in this readme enable the following:

1) **Metaclustering and Feature Engineering [Task 1]**: Run VoPo clustering on some FCS files and extract features (either frequency or functional)

2) **Using Engineered Features for Classification [Task 2]**: After extracting VoPo features, we show you how to apply feature selection and run the classification task

3) **Differentiation Score Visualization + Coloring by Phenotypic Marker Expression [Task 3]**: Generate comprehensive single-cell visualizations of frequency differences. We will show a generic example for general usage and examples for each of the 3 datasets in the paper (Figure 2A-C).  

4) **More Flexible: Visualize Cells According to a Differentiation Score Reflecting Functional Differences [Task 4]**: Examples on how to create single-cell visualizations with points colored by differentiation score (Fig 2 a-c) for frequency features and for functional features. We also have an option for 'directional differences' (e.g. to visualization which class a particular cell-type is higher in)

## Dependencies

Dependencies: We use a number of R packages. Please make sure that you have these packages installed. `flowCore`, `FNN`, `foreach`, `doParallel`, `iterators`, `plyr`, `randomForest`, `matrixStats`,`ROCR`,`FastKNN`, `miscTools`, `ggplot2`, `reshape2`, `viridis`, `pROC`, `igraph`.

## Installation Instructions

* FCS files in this tutorial are stored with `git-lfs`. Please make sure you have this installed (https://git-lfs.github.com/). 

* We used R version 3.4.4 to produce the results in this paper. R and the above packages are the only things that you need to install and we expect this to take minutes on a standard desktop computer.

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

Given FCS files, the first step is to create the metaclustering result. In this example, we the VoPo clustering object is named `Build` and can be generated as follows:

```R
#inputs in the order that they appear
	#S: The number of metaclustering iterations
	#K: The number of metaclusters per iteration
	#FileNames: The vector of filenames (with the whole path specified) e.g. FNames=c(~/FCS_1.fcs,~/FCS_1.fcs......)
	#doCPF='specify' means that we are going to specify the # of clusters per file
	#numCPF: The number of clusters per FCS file
	#MN: The vector of marker names that correspond to the columsn of the FCS files eg markerNames=c('CD45RA','CD3'.....)
	#ToUse: The vector of indices corresponding to the columns we want to use for clustering. For example, if you are clustering based only on phenotypic 		markers and they occur at indices 1, 3, and 5, we could define TU=c(1,3,5)
	#numCore: The number of cores to use for parallelization

Build=runRepMetaclust(S=50,K=50,FileNames=FNames,doCPF='specify',numCPF=1000,MN=markerNames,ToUse=TU,numCore=5)
```
Alternatively, please see a Demo Example (`Demo_VoPo.R`) for a collection of FCS files that will help to clarify the format of inputs

```R
source('Demo_Data/Demo_VoPo.R')
```

### Extracting Frequency Features

* If you want to extract frequency-related features from VoPo clustering object, `Build`, and the vector of filenames (called `FNames`) in the above `Demo_VoPo.R` script, you can do this with the following:

```R
source('VoPo_main/getFrequencyFeature.R')
FrequencyFeatures=getFrequencyFeature(Build,FNames)
```

This is the data matrix you can use for classification tasks (like we did in the paper). However, in the below example I show you how you can also get function-related features.

### Extracting Functional Features

* You can also extract functional features. You will input `Build` (the VoPo clustering object), `FNames`, which are the filenames (in the same order) as given to runRepMetaclust. Here we show results where the functional marker indices are given by `FInds`.

```R
source('VoPo_main/getFunctionalFeature')
#Here are the indices of the functional markers
FInds=readRDS('Processed/FI_HSR.rds')
FunctionalFeatures=getFunctionalFeature(Build,FNames,FInds)
```

**Note that you can use `FrequencyFeatures` or `FunctionalFeatures` either independently or concatenated with your favorite classifier.**

In task (2), we have an example with random forest, but you can feel free to use whatever you want with these features. We find that VoPo features + some regularized regression approach like Lasso or Elastic Net works very well for continuous outcomes.


## Task 2: Using Engineered Features for Classification.

We will show an example of how to use VoPo generated frequency features extracted using the `getFrequencyFeature.R` function. Our proposed example starts with a pre-processed VoPo clustering object obtained from the HSR dataset by running the function `runRepMetaclust.R`. 

Read in processed VoPo clustering result.

```R
Build_Surgery=readRDS('Processed/Build_Surgery.rds')
```

Load sample metadata, which we will use to extract the class labels for each sample and to keep samples from the same patient together during cross validation.

```R
Meta_Surgery=readRDS('Processed/Meta_Surgery.rds')
#The next line just gives the patient IDs for each sample
NameVec=Meta_Surgery$FileNames
```
Extract the frequency features

```R
FreqDF=getFrequencyFeature(Build_Surgery,NameVec)
```

We will generate a distribution of classification accuracies, repeating the cross validation 50 times. Note that the unsupervised feature selection approach also takes place here. Here is a brief description of the inputs

* FuncDF: Your feature matrix generated by `getFrequencyFeature.R`
* Y: vector of class labels (as factor) corresponding to the rows of FuncDF
* FPV: The number of features to use per clustering solution
* IterNumClus: Number of clusters per clustering iteration, extracted from the VoPo clustering object generated with `runRepMetaclust.R`
* propTrain: Proportion of patients to use for training. Here we are using 70%.
* numPerm: The number of leave-group-out cross validation trials to do 
* sampId: The IDs of the samples corresponding to the rows of FuncDF
* numCore: number of cores to use for parallelization

```R
ClAcc=runClassif(FuncDF=FreqDF,Y=Meta_Surgery$Class,FPV=40,IterNumClus=Build_Surgery$IterNumClus,propTrain=0.7,numPerm=50,sampID=as.character(Meta_Surgery$Subject),numCore=10)
```

## Task 3: Comprehensive Single-Cell Visualizations for each Dataset (Fig 2A-C.)

### General Usage

As described in the paper, the general idea is to sample a subset of cells across all sample FCS files, project them in 2D (for visualization purposes), and to map the differentiation scores onto them. 

We will show you the sequence of steps to use to apply this to your VoPo clustering result. In this example, we use tSNE for visualization, but you are free to construct your layout of cells however you want:

```R
#step 1: Sample cells across all FCS files

source('VoPo_main/SampleCells.R')

#Your inputs are 
	#FileNames: The vector of filenames that you used for runRepMetaclust.R
	#MN: the full vector of marker names (see Task 1)
	#ToUse: the indices of the markers that you used for clustering
	#NumCells: the number of cells to sample per file. We find 1000 works well.
	#NumCellsFinal: the number of cells to retain total (sampled randomly). We find 3000 works well

CellMat=SampleCells(FileNames,MN,ToUse,NumCells=1000,NumCellsFinal=30000)
```

Now, project cells obtained from the output of `SampleCells` using your favorite dimensionality reduction algorithm.

```R
library('Rtsne')
tRes=Rtsne(CellMat)$Y
```

### Color By The Expression of Each Phenotypic Marker

You can also use the following function to plot all of the cells in `CellMat` by the expression of each phenotypic marker for annotation.

```R
source('VoPo_main/vizAtlas_Phenotype')
```

### Coloring Cells by Differentiation Score
You are now ready to compute the per-cell differentiation score and to make visualizations

```R
source('VoPo_main/vizAtlas.R')

#Your inputs are
	#CellMat obtained above
	#Build is the VoPo clustering result obtained in Task (1)
	#Y is the vector of sample classes corresponding to the input files to VoPo clustering in task 1
	#ToUse is the indices of markers you used for clustering. These should be the same as your input to `SampleCells.R`
	#numCore is the number of cores to use
	#layout is the layout you computed on `CellMat` 
	#ourdir: the directory (`saveDir`) you will save the plot to. 

Atlas=vizAtlas(CellMat,Build,Y=myLabels,ToUse_Stroke,SampsToUse=NULL,numCore=35,layout=tRes,outdir=saveDir)
```
The above showed the sequence of steps you would need to compute a frequency-based differentiation score.

**Now, we show examples from the 3 datasets in the paper.** 

* You can look at any of these 3 examples to see how to use a VoPo clustering result for the comprehensive immune atlas visualization
* Results for each dataset will be within their respective folder in OutDir. 
* There is a plot for each surface marker showing its expression across cells (ex. CD3.jpg shows CD3 expression)
* pval.jpg plot colors the cells by their differentiation scores computed by the algorithm. 

### Normal Term Pregnancy Dataset (NTP)

```R
source('PaperFigures/Visualization/PregnancyViz.R')
```
You can now find plots for all markers and differentiation scores in OutDir/Pregnancy_Viz

### Hip Surgery Recovery Dataset (HSR)

Note in this example, we show how you can do visualization for a subset of samples by changing `SampsToUse=NULL`

```R
source('PaperFigures/Visualization/SurgeryViz.R')
```
You can now find plots for all markers and differentiation scores in OutDir/Surgery_Viz

### Longitudinal Stroke Recovery Dataset (LSR)

```R
source('PaperFigures/Visualization/StrokeViz.R')
```

You can now find plots for all markers and differentiation scores in OutDir/Stroke_Viz

## Task 4: Vizualize Differences Between Clinical Outcome Classes Based on Frequency or Function 

* Here we will show examples of how to color by differentation score for Frequency Based Features and Functional Features. **Note that in the paper we showed only examples for frequency-based features, but here we provide an example for function as well!**
* I will tell you the functions to use here. **More Specific Examples are coming soon!** Inputs are similar to what you have seen above.
* A differentiation score plot based on frequency will ultimately create a single plot, the same as that shown in Task 3.
* Using functional marker features, we will get one differentiation score plot for each functional marker.

Instructions for making differentiation score plots for both frequency and function

* Step 1: Select a large subset of cells across all of your FCS files. Use function, 

```R
source('VoPo_main/SampleCells.R') 
```

* Step 2: Color the sampled cells by phenotypic marker expression

```R
source('VoPo_main/vizAtlas_Phenotype.R')
```
* Step 3: Make a frequency map. Each cell will be colored by differentiation score according to frequency. Use the function,

```R
source('VoPo_main/vizAtlas_Freq.R')
```

* Step 4: Make functional maps. For each functional marker in your panel, each cell will be colored according to differentiation score based on the expression of that particular functional marker. Note that as one of the functional arguments you will specify the directory to write these plots to. They will all go to one place. Use the function,

```R
source('VoPo_main/vizAtlas_Function.R')
```

We also have functions for making directrional plots, which will indicate which class has higher frequency or function. By default, it will use the minimum value of the response vector (e.g. sample classifications) to be colored blue and the maximum value of the response vector (e.g. sample classifications) to be colored blue. So, if a point is colored red, it means it likely belongs to a cell-population where the functional marker expression/frequency was likely increased.

To make a frequency single-cell map with colors indicating direction:

```R
source('VoPo_main/vizAtlas_Freq_Directional.R')
```

To make a frequency single-cell map with colors indicating direction:

```R
source('VoPo_main/vizAtlas_Function_Directional.R')
```

This completes instructions for tasks 1-4 for general usage purposes.

__________________________________________________________________________

We also have examples (Tasks 5-6) for reproducting results in our paper.

5) **Task 5: Classification in all 3 datasets (Fig 2D)**: Feature Use VoPo features for classification tasks in each of the 3 datasets. Also reproduce the results in Figure 2D. Generating distributions of classification accuracies from single vs. repeated metaclustering solutions (Figure 2D.). 

Alternatively, in (6) we have provided an example for re-running all clustering results from scratch. However this will require downloading FCS files from flow repository and modifying the paths in the scripts to the data, accordingly.  

6) **Task 6: Re-run repeated metaclustering from scratch on the 3 datasets**: Re-running the repeated metaclustering strategy from scratch along with classification and visualization. Note that this task requires downloading FCS files from each dataset. 

-------------------------------------------------------------------------------------------------------------------------------

## Task 5: Classification Accuracies in all 3 Datasets (Reproduce Fig 2D.)

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

## Task 6: Re-run clustering from scratch. Use that clustering result to generate classification results and visualizations (Fig 2D.)

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



