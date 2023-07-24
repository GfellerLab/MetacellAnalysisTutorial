# Downstream analysis of metacells (for a discrete dataset)

Here we use the obtained metacell to run the downstream analysis on them instead of single-cell data. In this analysis, we treat metacell as single cell, neglecting information about their size (i.e., number of containing single cells). If you are interested in sample-weighted analysis, where metacell size is taken into account, see section \@ref(weighted-analysis).

## Standard analysis (R)
Standard analysis includes dimensionality reduction, clustering, differential expression etc using Seurat [ref] framework.



[//]: # (Standard downstream analysis of metacells using R (Seurat))

Under construction...

<!-- ### Dimensionality reduction  -->
<!-- ```{r supercell-2-seurat, eval=FALSE} -->
<!-- # MC.seurat <- supercell_2_Seurat( -->
<!-- #   SC.GE = MC.GE,  -->
<!-- #   SC = MC,  -->
<!-- #   fields = c("cell_line", "purity"), -->
<!-- #   var.genes = MC$genes.use, -->
<!-- #   N.comp = 10 -->
<!-- # ) -->

<!-- DimPlot(MC.seurat, cols = .color.cell.type, reduction = "umap") -->
<!-- ``` -->

<!-- ### Clustering  -->

<!-- ```{r r-mc-clustering, eval=FALSE} -->
<!-- MC.seurat <- FindClusters(MC.seurat, resolution = 0.5) -->
<!-- DimPlot(MC.seurat, reduction = "umap") -->
<!-- ``` -->

<!-- ### Differential expression analysis  -->

<!-- ```{r r-mc-DEG, eval=FALSE} -->
<!-- # Set idents to cell lines (as clusters are the same as cell lines) -->
<!-- Idents(MC.seurat) <- "cell_line" -->
<!-- levels(MC.seurat) <- sort(levels(Idents(MC.seurat))) -->

<!-- # Compute upregulated genes in each cell line (versus other cells) -->
<!-- MC.seurat.all.markers <-  FindAllMarkers( -->
<!--   MC.seurat,  -->
<!--   only.pos = TRUE, -->
<!--   min.pct = 0.25,  -->
<!--   logfc.threshold = 0.25,  -->
<!--   test.use = "t" -->
<!-- ) -->
<!-- saveRDS(MC.seurat.all.markers, file = file.path(data.folder, "output", paste0("MC_gamma_", gamma, "_all_markers_seurat.Rds"))) -->

<!-- # Top markers (select top markers of each cell line) -->
<!-- MC.seurat.top.markers <- MC.seurat.all.markers %>% -->
<!--    group_by(cluster) %>% -->
<!--     slice_max(n = 2, order_by = avg_log2FC) -->

<!-- MC.seurat.top.markers -->
<!-- ``` -->

<!-- ### Plot the expression of some markers -->

<!-- ```{r r-mc-clustering, eval=FALSE} -->
<!-- MC.seurat <- FindClusters(MC.seurat, resolution = 0.5) -->
<!-- DimPlot(MC.seurat, reduction = "umap") -->
<!-- ``` -->

<!-- ### Plot gene-gene correlation at single-cell and metacell levels  -->


<!-- ```{r r-mc-clustering, eval=FALSE} -->
<!-- MC.seurat <- FindClusters(MC.seurat, resolution = 0.5) -->
<!-- DimPlot(MC.seurat, reduction = "umap") -->
<!-- ``` -->


## Standard analysis (Python) {#standard-analysis-Py}


[//]: # (Standard downstream analysis of metacells using Py (Scanpy))

Standard analysis includes dimensionality reduction, clustering, differential expression etc using [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/#) framework.

Under construction...

### Dimensionality reduction 

### Clustering 

### Differential expression analysis 



## Advanced analysis

### GRN


[//]: # (GRN for metacells)

Gene correlation analysis suffers from large dropout rate of single-cell data and at the same time is very time and memory demanding. Metacells simultaneously adress both issues and thus are beneficial for gene co-expression and gene regulation analysis. Here we demonstrate usage of metacells for GRN analysis using SCENIC [ref]. 



## Sample-weighted analysis {#weighted-analysis}


[//]: # (Sample-weighted analysis of metacell using SuperCell framework)


One of the features of metacells are their size, which is a number of single cell it contains. Since metacells aggregate different number of cells, they also carry different amount of information. And thus, to better reproduce single-cell analysis, bigger metacells should have larger impact on the results than smaller metacells. For this, a sample-weighted analysis can be applied. Sample-weighted analysis is impleented in the SuperCell package. 




