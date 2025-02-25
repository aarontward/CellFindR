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

You should have your own datasets to run CellFindR but if you want to test it with our E14.5 Mouse Cochlea data it is accessible via geo: https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSM4037811

Here are the list of functions that pertains to CellFindR and the functions that helped generate figures from our paper. 

### load_tenx
This function loads 10x matrix files and goes through the general processing steps of Seurat and creates a usable Seurat object. Outputs a Seurat object

### out_tenx
This function is a quick function that selects a resolution, which then runs the UMAP and saves this as a pdf and as a Rdata Seurat object

### is_cluster
This function takes in a Seurat object and initialized cluster groupings and asks if the current clustering configuration satisfy's CellFindR's metrics. 
Output a boolean. 

### find_res
This function finds the highest resolution that satifies CellFindR metrics. Using an iterativ approach to the is_cluster function. This outputs a numeric value. 

### sub_clustering
This function takes in a Seurat object with initialized first layer clustering in the active ident and will output a Seurat object with CellFindR clustering groups. This function will take the most time to run. 
Once this is run, inside your directly should have a folder for each subcluster, such as 0, 1 and as well as folders for subclusters such as 1.0, 1.0.0 and etc. A RDS file with the cluster of interest. If there are subgroups within this group,  within each of these folders, it should include a matrix in the format of the get_matrix function as described below between the clusters, folders cluster and violin for the top differential expressed genes and a UMAP representation of this cluster. 

### get_analysis
For given Seurat object, returns mitochondria gene $ and UMI plots and the matrix, and statistics matrices in .csv form. Then also runs the get_matrix and get_stats functions below

### gen_matrix_plot
For given Seurat object, creates new directory inside the file and returns violin/cluster plots of the top 20 differential genes for the active ident. Will also return RDS Seurat file, generate a get_matrix matrix as below. 

### get_matrix
For Seurat object with clusters in active ident, will generate a matrix that shows average expression, differential gene expression and P-value associated with of the cluster for every gene. 
Across the rows, will be every gene name. Across the columns, the rows will be separated by the inputted identity cluster. If using the standard Seurat louvain clustering, it should start 0_Mean, 0_Avg_Diff, 0_Pval and so on for the other clusters. The mean reflects the log average that is centered by ScaleData output from the raw matrix file after Seurat standardization and calculated via the Seurat AverageExpression function. The Avg_diff computed via Wilcoxon rank sum test between the cluster of interest vs all the other cells not in this cluster. The Pval is the adjusted p-value, based on bonferroni correction using all genes in the dataset. These last two values are calculated from the FindMarker function in Seurat. 

### get_stats
For Seurat object with clusters in the active ident, will generate matrix that shows the cell count, nUMI and gene count for every cluster as well as the top expressed genes. 
Each column will be for one of the active identity in the object. The rows will show the metadata as the cell number, average number of reads, average number of numi in the cluster. THen it will show the top 10 differential genes as calculated via FindMarkers for this cluster compared to all the other cells not in this cluster in the object. 

## Running CellFindR vignette
Please see attached CellFindR_vignette.html or the RMarkdown file. This should walk you through the steps needed to generate CellFindR files. We will list out a list of files that are generated here at the end of the tutorial with the correct outputs.

At the end of the vignette, you should have an RDS file of the cellfindR run Seurat object. A folder that includes the files as run by get_analysis plot. The subclustering output of the folders. A cellfindR matrix (run through get_matrix). A folder with cluster and violin plots of top genes. 


#### Reference, Data Availability, and Contact:
If you find CellfindR useful, please reference:

Development of the Mouse and Human Cochlea at Single Cell Resolution. Kevin Shengyang Yu, Stacey M. Frumm, Jason S. Park, Katharine Lee, Daniel M. Wong, Lauren Byrnes, Sarah M. Knox, Julie B. Sneddon, Aaron D. Tward. bioRxiv 739680; doi: https://doi.org/10.1101/739680

If you have any inquiries or debugging questions, please contact Kevin Yu (yukev@uw.edu) or (aaron.tward@ucsf.edu)
