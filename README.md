# Advanced Exploratory Descriptive Analysis with Gene Expression Data
Small project on conducting advanced exploratory descriptive analysis (EDA) on a gene expression dataset.

## Problem

## Data
The data used from this projects come from a [this](https://archive.ics.uci.edu/ml/datasets/gene+expression+cancer+RNA-Seq) UCI - Machine Learning Repositry dataset with the following description reported in the website:  
  
This collection of data is part of the RNA-Seq (HiSeq) PANCAN data set, it is a random extraction of gene expressions of patients having different types of tumor: 
  
* **BRCA**: Breast Invasive Carcinoma 
* **KIRC**: Kidney Renal Clear Cell Carcinoma
* **COAD**: Colon Adenocarcinoma
* **LUAD**: Lung Adenocarcinoma
* **PRAD**: Prostate Adenocarcinoma

### Characteristics
Samples (instances) are stored row-wise. Variables (attributes) of each sample are RNA-Seq gene expression levels measured by illumina HiSeq platform. A dummy name (gene_XX) is given to each attribute, the attributes are in the same order as the [original submission](https://www.synapse.org/#!Synapse:syn4301332) where a complete list of the probes names can be found.  

**Number of instances**: 801  
**Number of attributes**: 20531  
**Missing Values**: None  
  
**Gene Expression Table**
|            | **gene_0** | ... | **gene_20530** |
|:----------:|:----------:|-----|----------------|
|  sample_0  |      0     | ... | 0.6            |
|     ...    |     0.6    | ... | 0              |
| sample_800 |     0.3    | ... | 0.2            |  
  
**Labels Table**  
|            | **Class** |
|:----------:|:---------:|
|  sample_0  |    PRAD   |
|     ...    |    ...    |
| sample_800 |    LUAD   |


### Source
Samuele Fiorini, University of Genoa, redistributed under [Creative Commons license](http://creativecommons.org/licenses/by/3.0/legalcode) from [here](https://www.synapse.org/#!Synapse:syn4301332)


## Analysis

### Pre-Processing

1. Retaining only the top 2000 genes showing the highest variance in the dataset.
2. Robust rescaling the expression levels of each gene, applying the formula `rescaled = (gene_expression - median(gene_expression)) / IQR(gene_expression)` where `IQR` stands for `Inter Quartile Range`.

### Principal Component Analysis (PCA)

<p align="center">   
  <img width="400" height="400"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/1.png">
</p>

<p align="center">   
  <img width="400" height="400"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/2.png">
</p>

### Uniform Manifold Approximation and Projection (UMAP) Analysis

<p align="center">   
  <img width="400" height="400"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/3.png">
</p>

### Density Based Clustering in UMAP Space Analysis

<p align="center">   
  <img width="400" height="400"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/4.png">
</p>

#### Sub-clusters Analysis

<p align="center">   
  <img width="700" height="200"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/13.png">
</p>

### Clusters Gene Expression Analysis 

<p align="center">   
  <img width="300" height="600"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/5.png">
</p>

<p align="center">   
  <img width="300" height="600"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/6.png">
</p>

### Lasso Regression Analysis 

<p align="center">   
  <img width="600" height="300"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/7.png">
</p>

<p align="center">   
  <img width="600" height="300"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/8.png">
</p>

<p align="center">   
  <img width="600" height="200"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/9.png">
</p>

<p align="center">   
  <img width="700" height="300"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/10.png">
</p>

<p align="center">   
  <img width="700" height="300"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/11.png">
</p>

### Cluster Comparison Analysis

<p align="center">   
  <img width="300" height="800"src="https://github.com/vb690/gene_expression_advanced_eda/blob/master/results/figures/12.png">
</p>
