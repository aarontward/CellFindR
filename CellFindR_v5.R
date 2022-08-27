# cellFindR 4.0, version for Seurat 4.0.5
#8.25.2022 Kevin Yu, Amar H Sheth

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(limma)
library(presto)

# loading tenx data
# input: file_loc = location of file to load
# res = resolution to run once
# proj_name = name of project
# cutoff = high end cut off of num genes
# output saves to file_loc a rdata file.
load_tenx <- function(file_loc, output_loc,proj_name = 'myproject',res = 0.1, cutoff= 10000, mito = FALSE){
  tenx.data <- Read10X(data.dir = file_loc)
  tenx <- CreateSeuratObject(counts = tenx.data, project = proj_name, min.cells = 3, min.features = 200)
  #processing
  tenx[["percent.MT"]] <- PercentageFeatureSet(tenx, pattern = "^MT-")
  VlnPlot(tenx, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  ###############################
  tenx <- NormalizeData(tenx, normalization.method = "LogNormalize", scale.factor = 10000)
  tenx <- FindVariableFeatures(tenx, selection.method = "vst", nfeatures = 2000,verbose=FALSE)
  top10 <- head(VariableFeatures(tenx), 10)
  plot1 <- VariableFeaturePlot(tenx)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  CombinePlots(plots = list(plot1, plot2))
  all.genes <- rownames(tenx)
  tenx <- ScaleData(tenx, features = all.genes, verbose = FALSE)
  #####################################
  tenx <- RunPCA(tenx, features = VariableFeatures(object = tenx), verbose = FALSE)
  tenx <- FindNeighbors(tenx, dims = 1:20,verbose=FALSE)
  return(tenx)
}

out_tenx <- function(tenx, output_loc,proj_name = 'myproject',res = 0.1, cutoff= 10000, mito = FALSE){
  tenx <- FindClusters(tenx, resolution = res,verbose=FALSE)
  tenx<-RunUMAP(tenx, dims = 1:20)
  ggsave(paste(output_loc, '/', proj_name, '_UMAP.pdf',sep = ""), 
         DimPlot(tenx, reduction = "umap"), width = 8, height = 8)
  #print("finished DIMPLOT")
  saveRDS(tenx, file = paste(output_loc,'/', proj_name, ".rds", sep = ''))
  return(tenx)
}

# asking if the grouping is a cluster
# input: tenx = tenx object
# thresh_genes = threshold of genes at thresh_val
# thresh_val = value of the threshold in log space
# pval = cut off of pval for the significance

is_cluster <- function(tenx, thresh_genes = 10, thresh_val = log(2), pval = 1e-4){
  val = 0 # groups that does not satisfy threshold genes
  counter = 0 # groups that satisfy threshold genes 
  # loop through the identitiy
  matrix_output <- data.frame(row.names = row.names(tenx))
  
  for (j in sort(unique(tenx@active.ident))){
    if (sum(tenx@active.ident == j) < 5){
      return(FALSE)
    }
    markers <- FindMarkers(tenx, ident.1 = j, min.pct = 0.25)
    markers <- markers[markers$p_val_adj < pval,]
    #markers <- wilcoxauc(tenx,"seurat_clusters","counts") %>% filter(group==j) %>% filter(pct_in>=25) %>%filter(padj<=pval_thresh)
    #find if the 10th biggest is less than log2, sum 
    #print(sort(markers$logFC, decreasing = TRUE)[thresh_genes])
    #find if the 10th biggest is less than log2, sum 
    #print(sort(markers$avg_log2FC, decreasing = TRUE)[thresh_genes])
    # if less than 10 significant genes
    if (length((markers$avg_log2FC)) < thresh_genes){
      val <- val + 1
    } else if (sort(markers$avg_log2FC, decreasing = TRUE)[thresh_genes] < thresh_val){
      #print(val)
      val <- val + 1
    } else{
      counter = counter + 1
    }
    if (val > 1){
      return(FALSE)
    }
  }
  if (val > 1){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

# finds resolution that satisfy
# input: tenx object
# initial resolution of starting clustering
# how much to increment up 
# threshold of genes
# value of the threshold 
find_res <- function(tenx, initial_res = 0.1, jump = 0.1, thresh_genes = 10, 
                     thresh_val = log(2)) {
  RES_POST <- initial_res #keeping
  RES_IT <- initial_res #iterative
  while(TRUE){
    tenx <- FindNeighbors(tenx, dims = 1:20,verbose=FALSE)
    tenx <- FindClusters(tenx, resolution = RES_IT,verbose=FALSE)
    
    # also check if theres only 1 cluster/ then can go up higher es
    # Find number of clusters
    length_group <- length(unique(tenx@active.ident))
    # if only one group then need to look deeper
    if (length_group == 1){
      # still not groups at 0.7 res stop and just step as 1
      if (RES_IT == 0.7){
        break
      }
    } else{
      testing <- is_cluster(tenx)
      if (testing == FALSE){ # if not real group
        #print(paste('go back', RES_IT, sep = ' '))
        RES_IT <- RES_IT - jump
        RES_POST <- RES_IT
        #print(RES_POST)
        break
      } else{ # valid groups
        RES_POST <- RES_IT
        #print(paste('ok',RES_IT, sep = ' '))
      }
    }
    RES_IT <- RES_IT + jump
  }
  # if there is only 1 group, return 0,
  return(RES_POST)
}

# getsubclustering
# input: tenx object
# location of output folder
# project_name
sub_clustering <- function(tenx, output_folder = '.', proj_name = 'proj_name',
                           thresh_genes = 10, thresh_val = log(2)){
  #print('Running subclustering')
  #public variables:
  ##resolution keeper:
  res_keep <- data.frame('cluster'= NA,'res'= NA, 'num_clusters' =NA)
  ##cellname cluster:
  celltocluster <- data.frame(row.names = colnames(tenx))
  celltocluster$cellnames <- colnames(tenx)
  
  #add cellfindr column to metadata as string
  tenx@meta.data$cellfindr <- as.character(tenx@active.ident)
  
  #what to iterate across:queue
  lib_c <- as.character(sort(unique(tenx@active.ident)))
  
  while(length(lib_c != 0)){
    # set to first value
    j <- lib_c[1]
    
    #print(paste('subsetting', j, sep = ''))
    sub_tenx <- subset(tenx, idents = toString(j))
    #print(paste('finished subsetting ', j, sep = ''))
    
    # need to recenter
    sub_tenx <- FindVariableFeatures(sub_tenx, selection.method = "vst", nfeatures = 2000,verbose=FALSE)
    sub_tenx <- FindNeighbors(sub_tenx, dims = 1:20, verbose = FALSE)
    if (dim(sub_tenx)[2] >49){
      sub_tenx <- RunUMAP(sub_tenx, dims = 1:20, n.neighbors = 10)
    }else{
      #print(paste("too small",j))
    }
    
    # get subgroups if there is a cluster, if not remove and label cells. 
    set_res <- find_res(sub_tenx)
    sub_tenx <-FindClusters(sub_tenx,pc.use = 1:20, resolution = set_res,verbose=FALSE)
    #print(paste("find clusters completed",j))
    # so subgroups, remove from queue: lib_c
    if (set_res == 0 || length(levels(sub_tenx@active.ident)) == 1){
      lib_c <- lib_c[lib_c != j]   
      #print(paste("if_statement_completed",j))
      file_create <-paste(output_folder,'/', j,'/',sep = '')
      dir.create(file_create)
      saveRDS(sub_tenx, file = paste(file_create,'/',j, '.rds', sep = ''))
    }
    # with subgroups
    else {
      # create the plots and matrix
      gen_matrix_plot(sub_tenx, output_folder, j)
      hold_group_names <- c()
      # find subgroups from j: 
      for (k in sort(unique(sub_tenx@active.ident))){
        l <- paste(j,k, sep = '.')
        # get column names to create
        sub2_tenx <- subset(sub_tenx, idents = toString(k))
        #print(paste("else_checkpoint",j))
        cellnames <- colnames(sub2_tenx)
        #print(paste("else_checkpoint2",j))
        rownames(tenx@meta.data) %in% cellnames
        #print(paste("else_checkpoint3",j))
        # add to the metadata file
        tenx@meta.data[rownames(tenx@meta.data) %in% cellnames,]$cellfindr = l
        #print(paste("else_checkpoint4",j))
        hold_group_names<- c(hold_group_names, l)
      }
      
      # add the new subgroups
      lib_c <- c(hold_group_names, lib_c)
      # remove original column
      lib_c <- lib_c[lib_c != j]   
    }
    #print(paste("else_statement_completed",j))
    tenx <- SetIdent(tenx, value = 'cellfindr')
  }
  
  #resort order of labels
  levels(tenx) <-str_sort(levels(tenx), numeric = TRUE)
  
  #graph umap
  ggsave(paste(output_folder, '/', proj_name, '_CellfindR_umap.pdf',sep = ""), 
         DimPlot(tenx, label = TRUE), width = 10, height = 8)
  return(tenx)
}

get_analysis <- function(tenx, output_folder = '.', proj_name = 'proj_name'){
  # output DataQ files
  
  # percent.mt plot
  ggsave(paste(output_folder, '/', proj_name, '_percent_mito.pdf',sep = ""), 
         FeaturePlot(tenx, "percent.MT"), width = 8, height = 8)
  
  # percent.mt violin plot
  ggsave(paste(output_folder, '/', proj_name, '_percent_mito_vln.pdf',sep = ""), 
         VlnPlot(tenx, "percent.MT"), width = length(levels(tenx@active.ident)), height = 5)
  
  # umi plot
  ggsave(paste(output_folder, '/', proj_name, '_nCount_RNA.pdf',sep = ""), 
         FeaturePlot(tenx, 'nCount_RNA'), width = 8, height = 8)
  
  # umi violin plot
  ggsave(paste(output_folder, '/', proj_name, '_nCount_RNA_vln.pdf',sep = ""), 
         VlnPlot(tenx, 'nCount_RNA'), width = length(levels(tenx@active.ident)), height = 5)
  
  #############
  #output Matrices
  
  tenx <-SetIdent(tenx, value = 'cellfindr')
  levels(tenx) <-str_sort(levels(tenx), numeric = TRUE)
  
  z <- get_matrix(tenx)
  write.csv(z, file = paste(output_folder, '/', 'matrix_cellfindr.csv', sep = ''))
  
  a <- get_stats(tenx)
  write.csv(a, file = paste(output_folder, '/', 'all_stats_cellfindr.csv', sep = ''))
  
  tenx <-SetIdent(tenx, value = 'seurat_clusters')
  y <- get_matrix(tenx)
  write.csv(y, file = paste(output_folder, '/', 'matrix_big_groups.csv', sep = ''))
  
  a <- get_stats(tenx)
  write.csv(a, file = paste(output_folder, '/', 'all_stats_big_groups.csv', sep = ''))
  
}


# generate matrix and plots
gen_matrix_plot <-function(tenx, output_folder = '.', proj_name = 'proj_name'){
  #create subdirectories
  file_create <- paste(output_folder, '/',proj_name, sep ='')
  dir.create(file_create)
  dir.create(paste(file_create, 'Cluster', sep = '/'))
  dir.create(paste(file_create, 'Violin', sep = '/'))
  ggsave(paste(file_create,'/',proj_name,'_umap.pdf', sep = ''), DimPlot(tenx, label = TRUE))
  #markers <- wilcoxauc(tenx,"seurat_clusters","counts") %>%  filter(logFC>=0.25) %>% filter(pct_in>=25)
  #markers_filtered <- markers %>% group_by(group) %>% top_n(n = 50, wt = logFC)
  #genes <- unique(markers_filtered$feature)
  markers <-FindAllMarkers(tenx,only.pos = TRUE,min.pct = 0.25,thresh.use = 0.25)
  markers_filtered <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
  genes <- unique(markers_filtered$gene)
  matrix_gen <- get_matrix(tenx)
  write.csv(matrix_gen,paste(file_create,'/',proj_name, '_matrix.csv', sep = ''), row.names = TRUE)
  saveRDS(tenx, file = paste(file_create,'/',proj_name, '.rds', sep = ''))
  
  for (i in genes) {
    # cluster maps
    ggsave(paste(paste(file_create, 'Cluster', '', sep = '/'),i,'.pdf',sep = ''),
           FeaturePlot(tenx, features = c(i),pt.size = 2), width = 8, height = 8)
    #violin plot
    ggsave(paste(paste(file_create, 'Violin', '', sep = '/'), i, '.pdf', sep = ''),    
           VlnPlot(tenx, c(i)),width = 6, height = 4)
  }
}

# generate matrix
get_matrix <- function(tenx){
  #print("getting matrix")
  avg_expression <- AverageExpression(tenx, slot = 'scale.data', verbose=FALSE) # log normalized, scale data
  matrix_all <- data.frame(row.names = rownames(avg_expression$RNA))
  for (i in levels(tenx@active.ident)){
    #print(i)
    #markers <- wilcoxauc(tenx,"seurat_clusters","counts") %>% filter(group==i) %>% filter(logFC>=0.1)
    #rownames(markers) = markers$feature
    markers <- FindMarkers(tenx, ident.1 = i, min.pct = 0.25)
    avg_val <- avg_expression$RNA[, toString(i)]
    avg_diff <- markers[rownames(avg_expression$RNA),]$avg_log2FC
    avg_diff[is.na(avg_diff)] <-0
    p_val <- markers[rownames(avg_expression$RNA),]$p_val_adj
    p_val[is.na(p_val)] <-1
    matrix_all <- cbind(matrix_all, avg_val)
    matrix_all <- cbind(matrix_all, avg_diff)
    matrix_all <- cbind(matrix_all, p_val)
  }
  name_col <- c()
  for (k in levels(tenx@active.ident)){
    #print(k)
    name_col <- c(name_col,(c(paste(k,'Mean', sep = '_'),paste(k,'Avg_diff', sep = '_') , paste(k,'Pval', sep = '_'))))
  }
  colnames(matrix_all) <- name_col
  return(matrix_all)
}
# get stats
get_stats <- function(tenx, num_genes = 20){
  aoe <- c("Group", "cell_number", "avg_read", "avg_umi")
  for (i in 1:num_genes){
    aoe <- c(aoe, paste('top_',i, sep = ""))
  }
  df <- data.frame(aoe)
  #initialize matrix
  for (groups in levels(tenx@active.ident)){
    subgroup <-subset(tenx, idents = groups)
    # group name
    aod <- c(groups)
    # cell number
    aod <- c(aod, length(subgroup@meta.data$nCount_RNA))
    # avg_read
    aod <- c(aod, mean(subgroup@meta.data$nCount_RNA))
    # avg_umi
    aod <- c(aod, mean(subgroup@meta.data$nFeature_RNA))
    # top 10 diff genes
    markers <- FindMarkers(tenx, groups)
    top_markers <- row.names(markers)[1:num_genes]
    for (topm in top_markers$feature){
      aod <- c(aod, topm)
    }
    df[groups] <-aod
  }
  return(df)
}

# getting plots
get_plots<- function(tenx, output_folder = '.'){
  dir_creater <- paste(output_folder, '/plots', sep = '')
  dir.create(dir_creater)
  for (groups in levels(tenx@active.ident)){
    markers <-FindMarkers(tenx, groups)
    maxer <- min(30, length(markers$feature))
    for (gene in markers$feature[1:maxer]){
      ggsave(paste(dir_creater,'/', gene, '_cluster.pdf', sep = ''),    
             FeaturePlot(tenx, features = gene), width = 6, height = 6)
      ggsave(paste(dir_creater,'/', gene, '_violin.pdf', sep = ''),    
             VlnPlot(tenx, features = gene, slot = "counts", log = TRUE), width = length(levels(tenx@active.ident)), height = 5)
    }
  }
}

# output metrics
metrics_output <- function(tenx, output_folder = '.', species = 'mouse'){
  tenx@active.assay
  # mito percent
  if (species == 'mouse'){
    ggsave(paste(output_folder, '/', 'percent_mito.pdf', sep = ''),    
           VlnPlot(tenx, features = 'percent.mt'), width =  length(levels(tenx@active.ident)), height = 5)
  }
  if (species == 'human'){
    ggsave(paste(output_folder, '/', 'percent_mito.pdf', sep = ''),    
           VlnPlot(tenx, features = 'percent.MT'), width =  length(levels(tenx@active.ident)), height = 5)
  }
  # umi
  #print('test')
  ggsave(paste(output_folder, '/', 'UMI.pdf', sep = ''),    
         VlnPlot(tenx, features = 'nFeature_RNA'), width =length(levels(tenx@active.ident)), height = 5)
}
