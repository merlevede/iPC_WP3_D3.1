# Deconvolution of multiple gene expression datasets in the framework of iPC project

We report here the code and approach used for the analyses performed in WP3 D3.1. 
The goal is to identify robust mechanims active in the tumor of interest in iPC, which are Ewing sarcoma, Medulloblastoma, Neuroblastoma, Hepatoblastoma and leukemias. We provide here the data needed to perform the analyses on Ewing Sarcoma tumor type. Due to size limitation, we provide only part of the datasets. Data for the other tumor types can be accessed through the following link (available as soon as they will be ready to use).

We advise the readers to start by reading the methods, then look at the example datasets and finally read the usage recommendations.

# Methods

The analysis of the collected datasets consisted in 4 steps, from decomposing a single dataset to interpreting signals found across many datasets. The workflow below summarizes the different steps and some details are given in the following sections. The approach applied here was already used in the study by Cantini et al. 2019: Assessing reproducibility of matrix factorization methods in independent transcriptomes.

Steps 1 and 2 were performed using an implementation in python of Stabilized ICA developped by Nicolas Captier, a member of Andrei Zinovyev team: https://github.com/ncaptier/Stabilized_ICA. We adapted this notebook to the needs of iPC project. 

## Step 1: sICA decomposition
The first step is based on matrix factorization, a classic way to reduce dimension of large datasets. There are 3 main dimensionality reduction methods: PCA, NMF and ICA. Based on our previous work on analyzing omics data, we selected sICA, an iterative application of ICA that ensures stability of the identified components. Each dataset is decomposed into 2 matrices, S containing the "metagenes" of dimension nb of genes * nb of components and A representing the "metasamples" of dimension nb of components * nb of samples.

We systematically applied sICA on each dataset using the following rule to determine the number of components, i.e. the dimension of A (in particular the number of columns) and S (in particular the number of rows!!! not with transpose): min(50, nb of samples/3). From our previous work, we show that usually 50 components is sufficient to see most of the relevant biological signals and taking too many components does not alter the quality of the first components (Kairov et al., 2017). 
Using this approach, we get a collection of metagenes and metasamples for each tumor type.

Step 1 was performed using python jupyter notebook sICA_notebook_ES_datasets. Inputs are normalized log data.

## Step 2: Reciprocal best hit analysis
Using different datasets of a same tumor type should allow to identify a same signal or pathway several times. Moreover, a process can be shared across different technologies or level of sequencing (e.g. bulk, single cell), or even between different data types (cell lines, PDX, patient or normal samples). Thus, an important and reliable mechanism should be identified, not in every, but in most of or several of the datasets. Less "contributive" processes are identified in less datasets.

To achieve this goal and identify robust communities, we took advantage of Reciprocal Best Hit (RBH) analysis. RBH was introduced by Cantini et al. 2019 in the context of omics data analysis. The idea is that there is at most one best match between any components of 2 datasets. Correlation is used between 2 identified metagenes being a best match to give a weight to this relationship. It indicates the confidence we can have in the fact that these 2 components coming from 2 independent datasets reflect the same mechanism.

Step 2 was performed using python jupyter notebook sICA_notebook_ES_datasets. Inputs are metagenes, generated at step 1.


## Step 3: Community detection and weighted metagenes
This step is divided into 2 parts. 

First, the detection of communities. Using results from RBH analysis, it is possible to identify groups reflecting the same processes across various datasets and collapse them into a community. A community is thus simply an ensemble of components sharing biological properties. 

To detect such communities, we used Markov Clustering algorithm (A. J. Enright et al., 2002 https://micans.org/mcl/). It automatically detects clusters formed by metagenes that are connected to each other. To consider a community as robust, i.e. a process found in several datasets, we required that it should contain at least 3 components. Of course, this parameter can by modified. 

The second part is the construction of the weighted metagenes, which are representatives of the communities by summarizing the information contained in their components. Knowing which metagenes are composing a community, we averaged their weights after making sure that the orientation of each metagene is correct. Indeed, in ICA, weights can be positive and negative. The weighted metagenes can now be furthered analyzed.

Step 3 was performed using R and the code in Rmd is provided here in mcl_and_weighted_metagenes.Rmd. Inputs are the results of reciprocal best hits, obtained at step 2. The json file should be converted to csv file, using cytoscape (https://cytoscape.org/) for example.


## Step 4: Community annotation
The ultimate goal is to characterize the identified robust communities. There are several ways to describe the communities:

1. Summarize the information at the level of a community: how many components does it contain? how many single cell / bulk datasets are there? how many patients / cell lines / PDX / controls datasets does it contain?
2. Identify the driver genes and perform over representation analyses on them, for both positive and negative tails.

We kept the genes beyond 3 S.D. from the mean of the weigths in each community.
We used clusterProfiler https://github.com/YuLab-SMU/clusterProfiler to performed over representation analyses using KEGG, Gene Ontology and Disease Ontology.

3. Perform gene sets enrichment analysis on the weighted metagenes

We used clusterProfiler https://github.com/YuLab-SMU/clusterProfiler for gsea

4. Perform metasample characterization when sample information are available

Step 4 was also performed using R and the code in Rmd is provided in Reports_communites_ES.Rmd and template_pathway_analysis.Rmd. Inputs are:

* the community composition obtained after the computation of weighted metagenes at step 3
* information about the datasets
* weighted metagenes, obtained at step 2, needed for gsea analyses

# Data

In this repository, we apply our approach to Ewing sarcoma, for which we gathered 20 datasets: 

* 5 bulk datasets of patients
* 4 single cell cell lines
* 8 single cell PDX
* 3 single cell controls (including 2 full embryos GSE137804)

These datasets log normalized are provided . Here we provide part of the datasets, due to github size limitation. The data are in folder data and zipped subfolder ES.

This is sufficient to test the method but we recommend to download the full dataset to get more insight into Ewing sarcoma tumor biologogy.


# Usage recommendations

## Step 1: sICA
To use the code provided here, please start by:

* unzipping the data in ES subdirectory. These are the inputs for step 1, that can be run through python notebook.
* indicating which data you would like to decompose by indicating the subfolder in the variable settings in the beginning of the script. 
* indicating what type of delimiter is used in the input gene expression files.

Currently, you can decompose all the files of one subfolder.

* providing the output path directory where you would like the results to be written

Outputs are:
* metagenes csv file
* metasamples csv file
* MDS pdf
* MSTD_coice pdf
* Stability_ICA_components pdf


## Step 2: RBH
Once you have run the matrix factorizations on the datasets you wanted, you can run the RBH analysis which is implemented in the last part of the python notebook.

Outputs are:
* RBH.json
* MNN pdf

## Step 3: Markov clustering and weighted metagenes construction
Step 2 generates a json file that should be converted to csv using cytoscape for example.

You can now run mcl_and_weighted_metagenes. The inputs needed are:

* rbh csv file
* metagenes files

Outputs are the weigthed metagenes. 

## Step 4: Community annotation 

In this last step, the annotation is performed by the description of each community.
Inputs needed are simply weigthed metagenes.

You can now run Reports_communites_ES.Rmd which uses template_pathway_analysis.Rmd. 
The output is an html report describing the communities.


