# Downstream analysis of metacells (for a continuous dataset)

Here we use the obtained metacell to run the downstream analysis on them instead of single-cell data. In this analysis, we treat metacell as single cell, neglecting information about their size (i.e., number of containing single cells). If you are interested in sample-weighted analysis, where metacell size is taken into account, see section \@ref(weighted-analysis-cont).

## Standard analysis (R)
Standard analysis includes dimensionality reduction, clustering, differential expression etc using Seurat [ref] framework.



[//]: # (Standard downstream analysis of metacells using R (Seurat))

### Dimensionality reduction 

### Clustering 

### Differential expression analysis 


## Standard analysis (Python) {#standard-analysis-Py}


[//]: # (Standard downstream analysis of metacells using Py (Scanpy))

Standard analysis includes dimensionality reduction, clustering, differential expression etc using [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/#) framework.

### Dimensionality reduction 

### Clustering 

### Differential expression analysis 





## Sample-weighted analysis {#weighted-analysis-cont}


[//]: # (Sample-weighted analysis of metacell using SuperCell framework)


One of the features of metacells are their size, which is a number of single cell it contains. Since metacells aggregate different number of cells, they also carry different amount of information. And thus, to better reproduce single-cell analysis, bigger metacells should have larger impact on the results than smaller metacells. For this, a sample-weighted analysis can be applied. Sample-weighted analysis is impleented in the SuperCell package. 




