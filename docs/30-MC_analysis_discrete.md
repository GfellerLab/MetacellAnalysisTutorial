# Downstream analysis of metacells 

In this chapter, we run standard and advanced downstream analyses on metacells instead of single-cell data. 
In this analysis, we treat each metacell as a single cell, neglecting information about the size of the metacell (i.e., number of containing single cells). 
If you are interested in sample-weighted analysis, where metacell size is taken into account, see section \@ref(weighted-analysis).

## Standard analysis (R) {#standard-analysis-R}



[//]: # (Standard downstream analysis of metacells using R (Seurat))



In this tutorial, standard analyses include dimensionality reduction, clustering and differential expression using the [Seurat](https://satijalab.org/seurat/) framework.


```r
library(Seurat) 
#> The legacy packages maptools, rgdal, and rgeos, underpinning this package
#> will retire shortly. Please refer to R-spatial evolution reports on
#> https://r-spatial.org/r/2023/05/15/evolution4.html for details.
#> This package is now running under evolution status 0
#> Attaching SeuratObject
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
```

### Load metacell Seurat object

We will use Seurat objects containing the metacells counts data and their annotation (*e.g.* cell-type annotation) and
proceed with standard Seurat downstream analyses.
Seurat objects containing metacells counts data and their annotation were generated at the end of sections \@ref(SuperCell-construction)
These objects can also be generated using the command line described in chapter \@ref(command-line)


```r
MC_tool = "SuperCell"
proj_name = "3k_pbmc"
annotation_column = "louvain"
celltype_colors <- c(
  "CD14+ Monocytes"    = "#E69F00",  # orange
  "B cells"            = "#56B4E9",  # sky blue
  "CD4 T cells"        = "#009E73",  # bluish green
  "NK cells"           = "#F0E442",  # yellow
  "CD8 T cells"        = "#0072B2",  # blue
  "FCGR3A+ Monocytes"  = "#D55E00",  # vermillion
  "Dendritic cells"    = "#CC79A7",  # reddish purple
  "Megakaryocytes"     = "#000000"   # black
)

MC.seurat = readRDS(paste0('./data/', proj_name, '/metacell_', MC_tool,'.rds'))
```


### Dimensionality reduction

As for single-cells, we normalize the raw counts (here aggregated raw counts) and we identify the most variable features in the metacells gene expression data.
Based on these features, we run PCA and use the first principal components to obtain a two dimensionnal representation of the data using UMAP.


```r
Idents(MC.seurat) <- annotation_column
MC.seurat <- NormalizeData(MC.seurat)
MC.seurat <- FindVariableFeatures(MC.seurat, selection.method = "vst", nfeatures = 2000)
MC.seurat <- ScaleData(MC.seurat)
#> Centering and scaling data matrix
MC.seurat <- RunPCA(MC.seurat, verbose = F)
MC.seurat <- RunUMAP(MC.seurat, dims = 1:30, verbose = F)
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session
DimPlot(MC.seurat, reduction = "umap", cols = celltype_colors)
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-dim-reduc-1.png" width="672" />

### Clustering

We cluster the metacells using Seurat clustering steps and visualize these clusters using UMAP:

```r
MC.seurat <- FindNeighbors(MC.seurat, dims = 1:30, reduction = "pca")
#> Computing nearest neighbor graph
#> Computing SNN
MC.seurat <- FindClusters(MC.seurat, resolution = 1.5)
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 264
#> Number of edges: 8154
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.5400
#> Number of communities: 6
#> Elapsed time: 0 seconds
DimPlot(MC.seurat, reduction = "umap", group.by = "seurat_clusters")
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-clustering-1.png" width="672" />

```r

cluster_colors=celltype_colors
names(cluster_colors)=c(0:7)
MC.seurat$SCclustering  <- SuperCell::supercell_cluster(D = dist(MC.seurat@reductions$pca@cell.embeddings[, 1:30]  ), k = 8)$clustering
DimPlot(MC.seurat, reduction = "umap", group.by = "SCclustering", cols = cluster_colors)
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-clustering-2.png" width="672" />

### Differential expression analysis

We perform diffrential analysis to identify the markers of CD8+ T cells as an example using the `FindMarkers` function. 

```r
# Set idents to metacell annotation 
Idents(MC.seurat) <- annotation_column
levels(MC.seurat) <- sort(levels(Idents(MC.seurat)))

CD8cells_markers <- FindMarkers(MC.seurat, ident.1 = "CD8 T cells", only.pos = TRUE)
#> For a more efficient implementation of the Wilcoxon Rank Sum Test,
#> (default method for FindMarkers) please install the limma package
#> --------------------------------------------
#> install.packages('BiocManager')
#> BiocManager::install('limma')
#> --------------------------------------------
#> After installation of limma, Seurat will automatically use the more 
#> efficient implementation (no further action necessary).
#> This message will be shown once per session
```

We visualize the top 5 markers for the CD8+ T cells.

```r
genes.to.plot <- rownames(CD8cells_markers)[1:5] 
VlnPlot(MC.seurat, features = genes.to.plot, ncol = 5, pt.size = 0.0, cols = celltype_colors)  
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-plot-genes-1.png" width="672" />

### Visualize gene-gene correlation 

We can use the `supercell_GeneGenePlot` function from the SuperCell package to visualize the correlation between marker genes of a cell-type: 
(i) at the single-cell level and
(ii) at the metacell level.

For that, we load the single-cell data from which the metacells were derived from.

```r
print(proj_name)
#> [1] "3k_pbmc"
sc_data <- readRDS(paste0("data/", proj_name, "/singlecell_seurat_filtered.rds"))
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize")
```

We visualize gene-gene correlation at the single-cell level:

```r

gene_x <- rownames(CD8cells_markers)[1:3]  
gene_y <- rownames(CD8cells_markers)[4:6]

alpha <- 0.7

p.sc <- SuperCell::supercell_GeneGenePlot(
  GetAssayData(sc_data, slot = "data"),
  gene_x = gene_x,
  gene_y = gene_y,
  clusters = sc_data@meta.data[, annotation_column],
  sort.by.corr = F, 
  alpha = alpha, 
  color.use = celltype_colors
)
p.sc$p
```

<img src="30-MC_analysis_discrete_files/figure-html/r-sc-plot-genes-1.png" width="672" />

We visualize gene-gene correlation at the metacell level:

```r

p.MC <- SuperCell::supercell_GeneGenePlot(GetAssayData(MC.seurat, slot = "data"),
                                          gene_x = gene_x,
                                          gene_y = gene_y,
                                          clusters = MC.seurat@meta.data[, annotation_column],
                                          sort.by.corr = F, supercell_size = MC.seurat$size,
                                          alpha = alpha,
                                          color.use = celltype_colors)
p.MC$p
```

<img src="30-MC_analysis_discrete_files/figure-html/r-MC-plot-genes-1.png" width="672" />


<!-- ## Standard analysis (Python) {#standard-analysis-Py} -->
<!-- ```{r, child='./sub_pages/30-downstream-Py.Rmd'} -->
<!-- ``` -->


## Sample-weighted analysis {#weighted-analysis}


[//]: # (Weighted downstream analysis of metacells using R (Seurat))



Code has been modified but text should be updated in this section!

```r
library(Seurat) 
library(dplyr)
library(ggplot2)
library(SuperCell)
```

### Load metacell Seurat object

We will use Seurat objects containing the metacells counts data and their annotation (*e.g.* and cell-type annotation) and
proceed with standard Seurat downstream analyses.
Seurat objects containing metacells counts data and their annotation were generated at the end of sections \@ref(SuperCell-construction)
These objects can also be generated using the command line described in chapter \@ref(command-line)


```r
MC_tool = "SuperCell"
proj_name = "3k_pbmc"
annotation_column = "louvain"
MC.seurat = readRDS(paste0('./data/', proj_name, '/metacell_', MC_tool,'.rds'))
```


### Dimensionality reduction

As for single-cells, we normalize the raw counts (here aggregated raw counts) and we identify the most variable features in the metacells gene expression data.
Based on these features, we run PCA and use the first principal components to obtain a two dimensionnal representation of the data using UMAP.


```r
MC.seurat <- NormalizeData(MC.seurat, normalization.method = "LogNormalize")

MC_list <- list(N.SC = ncol(MC.seurat),
                supercell_size = MC.seurat$size)
MC_list$PCA <- supercell_prcomp(
  Matrix::t(GetAssayData(MC.seurat, slot = "data")),
  genes.use = MC.seurat@misc$var_features,  # or a new set of HVG can be computed
  supercell_size = MC_list$supercell_size, # provide this parameter to run sample-weighted version of PCA,
  k = 10
)

MC_list$UMAP <- supercell_UMAP(
  SC = MC_list,
  PCA_name = "PCA",
  n_neighbors = 15 # large number to repel cells 
)

supercell_DimPlot(SC = MC_list, 
  groups = MC.seurat@meta.data[, annotation_column],
  dim.name = "UMAP", 
  title = paste0("UMAP of metacells colored by cell type assignment")
)
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-dim-reduc-weighted-1.png" width="672" />

### Clustering

We cluster the metacells using Seurat clustering steps and visualize these clusters using UMAP:

```r
# compute distance among metacells
D  <- dist(MC_list$PCA$x)

# cluster metacells
MC_list$clustering  <- supercell_cluster(D = D, k = 8, supercell_size = MC_list$supercell_size)

# Plot clustering result
supercell_DimPlot(
  MC_list, 
  groups = factor(MC_list$clustering$clustering),
  dim.name = "UMAP", 
  title = paste0("UMAP of metacells colored by metacell clustering")
)
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-clustering-weighted-1.png" width="672" />

### Differential expression analysis


```r
# Compute upregulated genes in each cell line (versus other cells)
MC.all.markers <- supercell_FindAllMarkers(
  ge = GetAssayData(MC.seurat, slot = "data"), 
  clusters = MC.seurat@meta.data[, annotation_column], 
  supercell_size = MC_list$supercell_size,
  only.pos = TRUE, 
  min.pct = 0, 
  logfc.threshold = 0.2
)

```

We select the top markers for each cell-type: 

```r
# Transform the output of `supercell_FindAllMarkers()` to be in the format of the `Seurat::FindAllMarkers()`
MC.all.markers.df <- data.frame()
for(cl in names(MC.all.markers)){
  cur <- MC.all.markers[[cl]]
  
  if(is.data.frame(cur)){
    cur$cluster <- cl
    cur$gene <- rownames(cur)
    cur$avg_log2FC <- cur$logFC
    MC.all.markers.df <- rbind(MC.all.markers.df, cur)
  } 
}

# Top markers (select top markers of each cell line)
MC.top.markers <- MC.all.markers.df %>%
   group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
```


We visualize the top 5 markers for the XX cells.

```r
Idents(MC.seurat) <- annotation_column
# genes.to.plot <- MC.seurat.top.markers$gene[MC.seurat.top.markers$cluster == unique(MC.seurat@meta.data[,annotation_column])[1]]
genes.to.plot <- MC.top.markers$gene[c(seq(1, 20, 5))]
VlnPlot(MC.seurat, features = genes.to.plot, ncol = 4, pt.size = 0.0)  
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-plot-genes-weighted-1.png" width="672" />




