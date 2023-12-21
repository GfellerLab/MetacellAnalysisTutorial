# Downstream analysis of metacells {#downstream-analysis}

In this chapter, we run standard and advanced downstream analyses on metacells instead of single-cell data. 
We will treat each metacell as a single cell, neglecting information about the size of the metacell (i.e., number of containing single cells). 
If you are interested in sample-weighted analysis, where metacell size is taken into account, see section \@ref(weighted-analysis).

## Standard analysis (R) {#standard-analysis-R}



[//]: # (Standard downstream analysis of metacells using R (Seurat))



In this tutorial, standard analyses include dimensionality reduction, clustering and differential expression using the [Seurat](https://satijalab.org/seurat/) framework.


```r
library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> The legacy packages maptools, rgdal, and rgeos, underpinning this package
#> will retire shortly. Please refer to R-spatial evolution reports on
#> https://r-spatial.org/r/2023/05/15/evolution4.html for details.
#> This package is now running under evolution status 0
#> 
#> Attaching package: 'SeuratObject'
#> The following object is masked from 'package:base':
#> 
#>     intersect
# If you have Seurat V5 installed, specify that you want to analyze Seurat V4 objects
wilcox.test <- "wilcox"
if(packageVersion("Seurat") >= 5) {
  options(Seurat.object.assay.version = "v4") 
  wilcox.test <- "wilcox_limma"
  print("you are using seurat v5 with assay option v4")}
#> [1] "you are using seurat v5 with assay option v4"
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
proj_name = "bmcite"
annotation_column = "celltype_simplified"

cell_types <- c("Prog_RBC", "Unconventional T", "Naive CD4 cell", "Non-Naive CD4 cell",
                "CD14 Mono", "B cell", "Naive CD8 cell", "Non-Naive CD8 cell",
                "NK", "GMP", "CD16 Mono", "pDC", "cDC2", "Prog_B 2",
                "Prog_Mk", "Plasmablast", "HSC", "LMPP", "Prog_DC", "Prog_B 1")

celltype_colors <- c("#7E57C2", "#1E88E5", "#FFC107", "#004D40", "#9E9D24",
                 "#F06292", "#546E7A", "#D4E157", "#76FF03", "#6D4C41",
                 "#26A69A", "#AB47BC", "#EC407A", "#D81B60", "#42A5F5",
                 "#2E7D32", "#FFA726", "#5E35B1", "#EF5350", "#3949AB")
names(celltype_colors) <-cell_types
MC.seurat = readRDS(paste0('./data/', proj_name, '/metacell_', MC_tool,'.rds'))
```


### Dimensionality reduction

As for single-cells, we normalize the raw counts (here aggregated raw counts) and we identify the most variable features in the metacells gene expression data.
Based on these features, we run PCA and use the first principal components to obtain a two dimensional representation of the data using UMAP.


```r
Idents(MC.seurat) <- annotation_column
MC.seurat <- NormalizeData(MC.seurat)
MC.seurat <- FindVariableFeatures(MC.seurat, selection.method = "vst", nfeatures = 2000)
MC.seurat <- ScaleData(MC.seurat)
#> Centering and scaling data matrix
MC.seurat <- RunPCA(MC.seurat, verbose = F)
MC.seurat <- RunUMAP(MC.seurat, dims = 1:30, verbose = F, min.dist = 1)
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session

data <- cbind(Embeddings(MC.seurat, reduction = "umap"),
              data.frame(size = MC.seurat$size,
                         cell_type = MC.seurat@meta.data[, annotation_column]))
colnames(data)[1:2] <- c("umap_1", "umap_2")
p_annot <- ggplot(data, aes(x= umap_1, y=umap_2, color = cell_type)) + geom_point(aes(size=size)) +
  ggplot2::scale_size_continuous(range = c(0.5,  0.5*max(log((data$size))))) +
   ggplot2::scale_color_manual(values = celltype_colors) +
  theme_classic() + guides(color=guide_legend(ncol=2))
p_annot
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-dim-reduc-1.png" width="672" />

```r
#DimPlot(MC.seurat, reduction = "umap", cols = celltype_colors, pt.size = log1p(MC.seurat$size))
```

### Clustering

We cluster the metacells using Seurat clustering steps and visualize these clusters using UMAP:

```r
MC.seurat$SCclustering  <- SuperCell::supercell_cluster(D = dist(MC.seurat@reductions$pca@cell.embeddings[, 1:30]  ), k = 15)$clustering
data <- cbind(Embeddings(MC.seurat, reduction = "umap"),
              data.frame(size = MC.seurat$size,
                         cluster = paste0("C_",MC.seurat$SCclustering)))
colnames(data)[1:2] <- c("umap_1", "umap_2")
p_cluster <- ggplot(data, aes(x= umap_1, y=umap_2, color = cluster)) + geom_point(aes(size=size)) +
  ggplot2::scale_size_continuous(range = c(0.5, 0.5*max(log1p((data$size))))) +
  theme_classic() + guides(color=guide_legend(ncol=2))
p_cluster
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-clustering-1.png" width="672" />

### Differential expression analysis

We perform differential analysis to identify the markers of our cluster 11 as an example using the `FindMarkers` function.

```r
# Set idents to metacell annotation
Idents(MC.seurat) <- "SCclustering"

cells_markers <- FindMarkers(MC.seurat, ident.1 = "11", only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, test.use = wilcox.test)
cells_markers[order(cells_markers$avg_log2FC, decreasing = T)[1:10], ]
#>                p_val avg_log2FC pct.1 pct.2    p_val_adj
#> SH2D1B  4.891817e-73   6.847440 0.976 0.063 8.320491e-69
#> CCNJL   7.548909e-60   6.408700 0.500 0.002 1.283994e-55
#> NCAM1   3.027264e-61   6.171451 0.833 0.052 5.149074e-57
#> KIR2DL1 5.624657e-40   5.907526 0.429 0.011 9.566979e-36
#> LGALS9C 1.299450e-31   5.633886 0.357 0.011 2.210234e-27
#> KLRF1   2.903180e-32   5.619850 1.000 0.406 4.938018e-28
#> NCR1    1.269200e-66   5.560546 0.952 0.069 2.158782e-62
#> LGALS9B 1.568102e-28   5.552238 0.286 0.006 2.667184e-24
#> SPON2   4.903245e-24   5.462028 0.929 0.469 8.339929e-20
#> NMUR1   7.401129e-57   5.457965 0.786 0.048 1.258858e-52
```

We see that the top marker genes for this cluster contain the killer cell lectin-like receptor (KLR) family,
which is a group of transmembrane proteins preferentially expressed in NK cells.


```r
genes = c("KLRF1", "KLRD1")
VlnPlot(MC.seurat, genes, ncol = 2, pt.size = 0.0)
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-plot-genes-1.png" width="672" />

We can verify the identification of the NK cell cluster by comparing the metacell annotation and the metacell clustering.


```r
p_cluster + p_annot
```

<img src="30-MC_analysis_discrete_files/figure-html/unnamed-chunk-5-1.png" width="960" />

### Visualize gene-gene correlation

We can use the `supercell_GeneGenePlot` function from the SuperCell package to visualize the correlation between marker genes of a cell-type:
(i) at the single-cell level and
(ii) at the metacell level.

For that, we load the single-cell data from which the metacells were derived from.

```r
print(proj_name)
#> [1] "bmcite"
sc_data <- readRDS(paste0("data/", proj_name, "/singlecell_seurat_filtered.rds"))
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize")
#> Warning: Command NormalizeData.RNA changing from SeuratCommand to SeuratCommand
```

We visualize gene-gene correlation at the single-cell level:

```r
cells_markers <- cells_markers[order(cells_markers$avg_log2FC, decreasing = T),]
gene_x <- rownames(cells_markers)[1:5]  
gene_y <- rownames(cells_markers)[6:10]

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
#> Warning: The `slot` argument of `GetAssayData()` is deprecated as of SeuratObject 5.0.0.
#> ℹ Please use the `layer` argument instead.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
p.sc$p
```

<img src="30-MC_analysis_discrete_files/figure-html/r-sc-gene_gene_cor-1.png" width="768" />

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

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-gene_gene_cors-1.png" width="672" />


## Standard analysis (Python) {#standard-analysis-Py}


[//]: # (Standard downstream analysis of metacells using Py (Scanpy))



In this section, standard analysis includes dimensionality reduction, clustering, differential expression etc using [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/#) framework.


```python
import os
import numpy as np
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 1 
```



```python
MC_tool = "SEACells"
proj_name = "bmcite"
annotation_column = "celltype_simplified"
adata = sc.read(os.path.join('./data', proj_name, f'metacell_{MC_tool}.h5ad'))
```

### Dimensionality reduction

```python
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Freeze the state of the AnnData object for later use in differential testing and visualizations of gene expression
adata.raw = adata # step needd only if I use regress_out steps

# Compute PCA (highly variable genes will be used)
sc.tl.pca(adata, svd_solver='arpack')

# Compute the neighbor graph
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
#> /opt/conda/envs/MetacellAnalysisToolkit/lib/python3.9/site-packages/umap/distances.py:1063: NumbaDeprecationWarning: [1mThe 'nopython' keyword argument was not supplied to the 'numba.jit' decorator. The implicit default value for this argument is currently False, but it will be changed to True in Numba 0.59.0. See https://numba.readthedocs.io/en/stable/reference/deprecation.html#deprecation-of-object-mode-fall-back-behaviour-when-using-jit for details.[0m
#>   @numba.jit()
#> /opt/conda/envs/MetacellAnalysisToolkit/lib/python3.9/site-packages/umap/distances.py:1071: NumbaDeprecationWarning: [1mThe 'nopython' keyword argument was not supplied to the 'numba.jit' decorator. The implicit default value for this argument is currently False, but it will be changed to True in Numba 0.59.0. See https://numba.readthedocs.io/en/stable/reference/deprecation.html#deprecation-of-object-mode-fall-back-behaviour-when-using-jit for details.[0m
#>   @numba.jit()
#> /opt/conda/envs/MetacellAnalysisToolkit/lib/python3.9/site-packages/umap/distances.py:1086: NumbaDeprecationWarning: [1mThe 'nopython' keyword argument was not supplied to the 'numba.jit' decorator. The implicit default value for this argument is currently False, but it will be changed to True in Numba 0.59.0. See https://numba.readthedocs.io/en/stable/reference/deprecation.html#deprecation-of-object-mode-fall-back-behaviour-when-using-jit for details.[0m
#>   @numba.jit()
#> /opt/conda/envs/MetacellAnalysisToolkit/lib/python3.9/site-packages/umap/umap_.py:660: NumbaDeprecationWarning: [1mThe 'nopython' keyword argument was not supplied to the 'numba.jit' decorator. The implicit default value for this argument is currently False, but it will be changed to True in Numba 0.59.0. See https://numba.readthedocs.io/en/stable/reference/deprecation.html#deprecation-of-object-mode-fall-back-behaviour-when-using-jit for details.[0m
#>   @numba.jit()

# Run umap
sc.tl.umap(adata)

# Plot metacells in the UMAP space
sc.pl.umap(adata, color=[annotation_column], size = 100)
#> /opt/conda/envs/MetacellAnalysisToolkit/lib/python3.9/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
#>   cax = scatter(
```

<img src="30-MC_analysis_discrete_files/figure-html/py-mc-dim-reduc-1.png" width="672" />

### Clustering

Perform clustering on the metacell data.

```python
# run laiden graph-clustering
sc.tl.leiden(adata, neighbors_key = "neighbors", resolution = 2)

# plot the metacells in the UMAP space and color by cluster
sc.pl.umap(adata, color=['leiden'], size = 100)
#> /opt/conda/envs/MetacellAnalysisToolkit/lib/python3.9/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
#>   cax = scatter(
```

<img src="30-MC_analysis_discrete_files/figure-html/py-mc-clustering-3.png" width="672" />

### Differential expression analysis

Identify marker genes for each group of metacells. 
We visualize top markers for NK Tcells (cluster 6).

```python
# Identify marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
#> /opt/conda/envs/MetacellAnalysisToolkit/lib/python3.9/site-packages/scanpy/plotting/_tools/__init__.py:397: UserWarning: Attempting to set identical low and high ylims makes transformation singular; automatically expanding.
#>   ax.set_ylim(ymin, ymax)
```

<img src="30-MC_analysis_discrete_files/figure-html/py-mc-markers-5.png" width="2688" />

```python

# Show the top marker genes
print(pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5))
#>          0        1        2      3  ...        7               8         9      10
#> 0    TPD52  PIK3IP1   S100A8   LAG3  ...     PKIB           TRADD     PILRA    SPIB
#> 1    CD79A     LEF1     MNDA  TRGC2  ...    NDRG2           CLIC5    STXBP2    CLN8
#> 2  POU2AF1    SARAF  C14orf2   LYAR  ...     ENHO  RP11-1399P15.1    CDKN1C   PTPRS
#> 3     LCN8     CCR7    APLP2  DUSP2  ...    GPAT3            AQP3  C19orf38  BCL11A
#> 4  RASGRP3  FAM134B    MARC1   IL32  ...  CLEC10A         TNFRSF4     COTL1    IRF8
#> 
#> [5 rows x 11 columns]

# Visualize marker genes
sc.pl.violin(adata, ['KLRF1', 'IL2RB', 'GNLY'], groupby=annotation_column, size = 2, rotation = 90)
```

<img src="30-MC_analysis_discrete_files/figure-html/py-mc-markers-6.png" width="2564" />


## Sample-weighted analysis {#weighted-analysis}


[//]: # (Weighted downstream analysis of metacells using R (Seurat))




```r
library(Seurat)
# If you have Seurat V5 installed, specify that you want to analyze Seurat V4 objects
if(packageVersion("Seurat") >= 5) {options(Seurat.object.assay.version = "v4"); print("you are using seurat v5 with assay option v4")}
#> [1] "you are using seurat v5 with assay option v4"
library(dplyr)
library(ggplot2)
library(SuperCell)
```

### Load metacell Seurat object

We will use Seurat objects containing the metacells counts data and their annotation (*e.g.* and cell-type annotation) and
proceed with downstream analyses considering the size of each metacells.
Seurat objects containing metacells counts data and their annotation were generated at the end of sections \@ref(SuperCell-construction)
These objects can also be generated using the command line described in chapter \@ref(command-line)


```r
MC_tool = "SuperCell"
proj_name = "bmcite"
annotation_column = "celltype_simplified"

cell_types <- c("Prog_RBC", "Unconventional T", "Naive CD4 cell", "Non-Naive CD4 cell",
                "CD14 Mono", "B cell", "Naive CD8 cell", "Non-Naive CD8 cell",
                "NK", "GMP", "CD16 Mono", "pDC", "cDC2", "Prog_B 2",
                "Prog_Mk", "Plasmablast", "HSC", "LMPP", "Prog_DC", "Prog_B 1")

celltype_colors <- c("#7E57C2", "#1E88E5", "#FFC107", "#004D40", "#9E9D24",
                 "#F06292", "#546E7A", "#D4E157", "#76FF03", "#6D4C41",
                 "#26A69A", "#AB47BC", "#EC407A", "#D81B60", "#42A5F5",
                 "#2E7D32", "#FFA726", "#5E35B1", "#EF5350", "#3949AB")
names(celltype_colors) <-cell_types
MC.seurat = readRDS(paste0('./data/', proj_name, '/metacell_', MC_tool,'.rds'))
```


### Dimensionality reduction

As for single-cells, we normalize the raw counts (here aggregated raw counts) and we identify the most variable features in the metacells gene expression data.
Based on these features, we run a sample weighted PCA using the function `supercell_prcomp` from the SuperCell R package
and use the first principal components to obtain a two dimensionnal representation of the data using UMAP.
Using the `supercell_DimPlot` function from the the SuperCell R package we can visualize the metacells and their sized in UMAP space.


```r
MC.seurat <- NormalizeData(MC.seurat, normalization.method = "LogNormalize")

MC_list <- list(N.SC = ncol(MC.seurat),
                supercell_size = MC.seurat$size)
MC_list$PCA <- SuperCell::supercell_prcomp(
  Matrix::t(GetAssayData(MC.seurat, slot = "data")),
  genes.use = MC.seurat@misc$var_features,  # or a new set of HVG can be computed
  supercell_size = MC_list$supercell_size, # provide this parameter to run sample-weighted version of PCA,
  k = 30
)

MC_list$UMAP <- supercell_UMAP(
  SC = MC_list,
  PCA_name = "PCA",
  n.comp = 30, n_neighbors = 15, min_dist=0.5
)

supercell_DimPlot(SC = MC_list,
  groups = MC.seurat@meta.data[, annotation_column],
  dim.name = "UMAP",
  title = paste0("UMAP of metacells colored by cell type assignment"), color.use = celltype_colors
) + guides(color=guide_legend(ncol=2))
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-dim-reduc-weighted-1.png" width="672" /><img src="30-MC_analysis_discrete_files/figure-html/r-mc-dim-reduc-weighted-2.png" width="672" />

### Clustering

We cluster the metacells using the function `supercell_cluster` from SuperCell R package to perform the clustering step and visualize these clusters in the UMAP space:

```r
# compute distance among metacells
D  <- dist(MC_list$PCA$x)

# cluster metacells
MC_list$SCclustering  <- supercell_cluster(D = D, k = 15, supercell_size = MC_list$supercell_size)
MC.seurat$SCclustering <- MC_list$SCclustering$clustering

# Plot clustering result
supercell_DimPlot(
  MC_list,
  groups = factor(MC_list$SCclustering$clustering),
  dim.name = "UMAP",
  title = paste0("UMAP of metacells colored by metacell clustering")
)  + guides(color=guide_legend(ncol=2))
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-clustering-weighted-1.png" width="672" /><img src="30-MC_analysis_discrete_files/figure-html/r-mc-clustering-weighted-2.png" width="672" />

### Differential expression analysis

We perform diffrential analysis to identify the markers of our clusters using the `supercell_FindAllMarkers` function from the SuperCell package.

```r
# Compute upregulated genes in each cell line (versus other cells)
MC.all.markers <- supercell_FindAllMarkers(
  ge = GetAssayData(MC.seurat, slot = "data"),
  clusters = MC_list$SCclustering$clustering,
  supercell_size = MC_list$supercell_size,
  only.pos = TRUE,
  min.pct = 0,
  logfc.threshold = 0.2
)
```

We select the markers for cluster 9:

```r
cluster_markers <- MC.all.markers[[9]]
MC.top.markers <- cluster_markers[order(cluster_markers$logFC, decreasing = T),]
head(MC.top.markers)
#>       p.value adj.p.value     pct.1     pct.2    logFC w.mean.1   w.mean.2
#> GNLY        0           0 1.0000000 0.9692845 3.444786 4.629992 0.54215411
#> NKG7        0           0 1.0000000 0.9536207 2.853815 4.023960 0.51434288
#> GZMB        0           0 0.9916981 0.6401282 2.529441 2.771860 0.15034752
#> KLRF1       0           0 1.0000000 0.5681106 2.225002 2.353682 0.09508412
#> KLRD1       0           0 1.0000000 0.5480537 2.184094 2.508815 0.17047908
#> CST7        0           0 1.0000000 0.7280256 2.081451 2.706358 0.30405936
```

We visualize the top 5 markers for the cluster 9 and see that the top marker genes for this cluster contain marker genes of natural killer cells such as GZMB and GNLY.

```r
Idents(MC.seurat) <- "SCclustering"
# genes.to.plot <- MC.seurat.top.markers$gene[MC.seurat.top.markers$cluster == unique(MC.seurat@meta.data[,annotation_column])[1]]
# genes.to.plot <- MC.top.markers$gene[c(seq(1, 20, 5))]
genes.to.plot <- rownames(MC.top.markers)[1:5]
VlnPlot(MC.seurat, features = genes.to.plot, ncol = 5, pt.size = 0.0)  
#> Warning: Groups with fewer than two data points have been dropped.
#> Groups with fewer than two data points have been dropped.
#> Groups with fewer than two data points have been dropped.
#> Groups with fewer than two data points have been dropped.
#> Groups with fewer than two data points have been dropped.
```

<img src="30-MC_analysis_discrete_files/figure-html/r-mc-plot-genes-weighted-1.png" width="960" />


