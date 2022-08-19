# CellFindR
scRNAseq tools for hierarchical cell clustering

## Description
CellFindR is an algorithm that iteratively applies Louvain clustering to single cell datasets in an automated fashion to identify biologically meaningful cell subpopulations.  It utilizes many of the functions within Seurat and is implemented as a wraparound on that package. CellFindR's outputs include matricies with gene expression information for each subpopulation as well as automated generation of violin plots and feature plots for the top differentially expressed genes.  It also includes additional data visualization functionalities.  CellfindR was developed by Kevin Yu and updated by Amar Sheth in the Tward lab at UCSF.

## Step 1: Downloading pre-requisites
CellFindR v4.0 works and has been tested in R version 4.0.5 and Seurat version 4.0.5. To download R and Seurat, please see links below:

- https://satijalab.org/seurat/articles/install.html
- https://cran.r-project.org/src/base/R-4/

Furthermore, CellFindR requires the following packages to be downloaded: dplyr, stringr, ggplot2, limma, and presto. To download dplyr, stringr, and ggplot2, you can simply type install.packages("package+name") in the R console. To download limma and presto see code chunks below.

Limma: https://kasperdanielhansen.github.io/genbioconductor/html/limma.html
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
```
Presto: https://github.com/immunogenomics/presto
```R
library(devtools)
install_github('immunogenomics/presto')
```

Downloading RColorBrewer and pheatmap is optional.

## Step 2: Downloading CellFindR
Clone the repository from this github: https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository. 

## Step 3: How to use
Using CellFindR involves first loading in the functions provided in the CellFindR.R script. The best way do this is to open the file from the repository folder, highlight all (Ctrl/Cmd + A), and Run (Shift + Enter). Once loaded, we suggest creating a new R script (File > New > R Script) and running commands listed and demonstrated in the CellFindR_vignette.html file. The vignette can be opened using any web browser. **We detail each of the functions below and suggestions for ways in which our code can be modified for your use.** Please note that you will have to rerun the modified function to update it before calling it. 

### load_tenx

### out_tenx

### is_cluster

### find_res

### sub_clustering

### get_analysis

### gen_matrix_plot

### get_matrix

### get_stats

### get_plots

### metrics_output

#### Reference, Data Availability, and Contact:
If you find CellfindR useful, please reference:

Development of the Mouse and Human Cochlea at Single Cell Resolution. Kevin Shengyang Yu, Stacey M. Frumm, Jason S. Park, Katharine Lee, Daniel M. Wong, Lauren Byrnes, Sarah M. Knox, Julie B. Sneddon, Aaron D. Tward. bioRxiv 739680; doi: https://doi.org/10.1101/739680

If you have any inquiries or debugging questions, please contact Amar Sheth (amar.sheth@yale.edu) or (aaron.tward@ucsf.edu)
