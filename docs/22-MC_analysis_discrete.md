# Downstream analysis of metacells (for a discrete dataset)

Here we use the obtained metacell to run the downstream analysis on them instead of single-cell data. In this analysis, we treat metacell as single cell, neglecting information about their size (i.e., number of containing single cells). If you are interested in sample-weighted analysis, where metacell size is taken into account, see section \@ref(weighted-analysis).

## Standard analysis (R)
Standard analysis includes dimensionality reduction, clustering, differential expression etc using Seurat [ref] framework.



[//]: # (Standard downstream analysis of metacells using R (Seurat))

Under construction...

### Dimensionality reduction 

### Clustering 

### Differential expression analysis 


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




