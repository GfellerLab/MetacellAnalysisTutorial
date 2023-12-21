# Requirements{.unnumbered}



This chapter describes how to obtain the packages and data needed to reproduce the analyses performed in this tutorial.

## Installations {#installations}

### Using conda (recommended)
To build a conda environment containing the three metacell building tools used in this tutorial (SuperCell, MC2 and SEACells), 
please follow the instructions provided in the README of our MetacellAnalysisToolkit [github repository](https://github.com/GfellerLab/MetacellToolkit).

Then run the following lines to define the python path to use.

```r
library(reticulate)
conda_env <-  conda_list()[reticulate::conda_list()$name == "MetacellAnalysisToolkit","python"]

use_condaenv(conda_env)
```

The following R packages should also be installed. 
This tutorial was developed under Seurat V4 but is also compatible with Seurat V5.


```r
remotes::install_github("GfellerLab/SuperCell",upgrade = "never")
remotes::install_github("GfellerLab/MetacellAnalysisToolkit",upgrade = "never")
remotes::install_github("rstudio/reticulate",upgrade = "never")  #temporary fix for reading sparse matrix with R anndata https://github.com/rstudio/reticulate/issues/141
#install.packages("Seurat") # uncoment to update Seurat to V5 since V5 not yet on conda
#BiocManager::install('limma',update = F) # uncoment if Seurat V5 used 
```

### Without conda
If you don't have conda, you can use the following instructions:

Set up a python virtual environment with MC2 and SEACells installed:


```bash
pip install virtualenv
virtualenv my_env
source my_env/bin/activate

# Installing SEACells
pip install git+https://github.com/dpeerlab/SEACells

# Install MC2
pip install git+https://github.com/tanaylab/metacells
```

In R, install the SuperCell package:

```r
remotes::install_github("GfellerLab/SuperCell", upgrade = "never")
```

To run python function in R, install reticulate:

```r
install.packages('reticulate')
```

To use the python libraries installed in the virtual environment, define the RETICULATE_PYTHON variable as follow:

```bash
echo 'RETICULATE_PYTHON=my_env/bin/python' > '.Renviron'
```

The following R packages should also be installed. 
This tutorial was developed under Seurat V4 but is also compatible with Seurat V5.


```r
remotes::install_github("GfellerLab/SuperCell",upgrade = "never")
remotes::install_github("GfellerLab/MetacellAnalysisToolkit",upgrade = "never")
remotes::install_github("rstudio/reticulate",upgrade = "never")  #temporary fix for reading sparse matrix with R anndata https://github.com/rstudio/reticulate/issues/141
#install.packages("Seurat") # uncoment to update Seurat to V5 since V5 not yet on conda
#BiocManager::install('limma',update = F) # uncoment if Seurat V5 used 
```

## Retrieve a discrete dataset (Bone marrow dataset) {#bmcite-data}

To test metacell construction on a discrete dataset, we retrieved the "bmcite" dataset from the SeauratData R package containing around 30'000 cells.

The data are saved in the following file for future analyses in R (use of SuperCell): "data/bmcite/singlecell_seurat_filtered.rds".


```r
library(SeuratData)
InstallData("bmcite")

data("bmcite")
bmcite
head(bmcite@meta.data)
bmcite$celltype_simplified <- plyr::revalue(bmcite$celltype.l2, 
                                            c("CD8 Effector_1" = "Non-Naive CD8 cell",
                                              "CD8 Effector_2" = "Non-Naive CD8 cell",
                                              "CD8 Memory_1" = "Non-Naive CD8 cell",
                                              "CD8 Memory_2" = "Non-Naive CD8 cell",
                                              "CD8 Naive" = "Naive CD8 cell",
                                              "CD4 Naive" = "Naive CD4 cell",
                                              "CD4 Memory" = "Non-Naive CD4 cell",
                                              "Treg" = "Non-Naive CD4 cell",
                                              "Naive B" = "B cell",
                                              "Memory B" = "B cell",
                                              "CD56 bright NK" = "NK",
                                              "MAIT" = "Unconventional T",
                                              "gdT" = "Unconventional T"
                                              ))
bmcite <- bmcite[,-grep("Prog",bmcite$celltype_simplified)]
if(packageVersion("Seurat") >= 5) {
  bmcite[["RNA"]] <- as(object = bmcite[["RNA"]], Class = "Assay")
}
saveRDS(bmcite, file = paste0("data/bmcite/singlecell_seurat_filtered.rds"))

```


The data are saved in the following file for future analyses in python (use of SEACells and MC2): "data/bmcite/singlecell_anndata_filtered.h5ad".

```r
library(anndata)
adata <- AnnData(X = Matrix::t(bmcite@assays$RNA@counts),
                 var = data.frame(genes = rownames(bmcite@assays$RNA@counts)),
                 obs = bmcite@meta.data)

write_h5ad(adata, paste0("data/bmcite/singlecell_anndata_filtered.h5ad"))

```

## Retrieve a continuous dataset (CD34 dataset) {#CD34-data}

To test metacell construction on continuous dataset, we retrieved the CD34 dataset provided in [@SEACells]:

```bash
mkdir -p data/CD34
wget -O data/CD34/cd34_multiome_rna.h5ad 'https://dp-lab-data-public.s3.amazonaws.com/SEACells-multiome/cd34_multiome_rna.h5ad' 
```

The downloaded file will be used in the section \@ref(command-line).

## Retrieve the lung atlas dataset {#HLCA-data}

This dataset will be used for the integration of a large number of single-cell datasets at the level of metacells (see section \@ref(integration)).
Considering, the large size of the data to download, if you don't consider running the integration analysis, you can skip this part of the tutorial.

### Downloading the atlas {#HLCA-data-download}

To illustrate how metacells can be used in the context of single-cell data integration,
we used a cell atlas of the human lung (core) available on [cellxgene](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293). 
To download the data, please choose the `.h5ad` option after clicking on the download button for the core atlas (3 tissues, 584'944 cells).

Save these data in the `data/HLCA/` directory. 

Please note that this may take some time (\~45 mins) as the file is quite large (5.6 GB).

###  Splitting atlas by datasets

We will use anndata to read in backed mode (saving a lot of memory) the whole atlas and write one h5ad file for each dataset. 
This should take less than 10 minutes.

If you are limited in time feel free to process only a subset of the dataset.


```r
t0.split <- Sys.time()

library(anndata)
adata <- read_h5ad("data/HLCA/local.h5ad",backed = "r")
adata$var_names <- adata$var$feature_name # We will use gene short name for downstream analyses
datasets <- unique(adata$obs$dat)

# If you are limited in time you can process on half of the datasets (uncomment th following line)
# datasets <- datasets[1:7]

print(dim(adata))

lapply(datasets,FUN =  function(x) {
  dir.create(paste0("data/HLCA/datasets/",x),recursive = T)
  adata.dataset <- AnnData(X = adata[adata$obs$dataset == x]$raw$X,
                           var = adata[adata$obs$dataset == x]$var,
                           obs = adata[adata$obs$dataset == x]$obs)
  #This will allow us to construct supervised metacell for each cell type in each sample later in the tutorial
  adata.dataset$obs$ann <- as.character(adata.dataset$obs$ann_level_3)
  # For cell without an annotation at the 3rd level we will use the second level of annotation
  adata.dataset$obs$ann[adata.dataset$obs$ann_level_3 == 'None'] = as.character(adata.dataset$obs$ann_level_2[adata.dataset$obs$ann_level_3 == 'None'])
  adata.dataset$obs$ann_sample <- paste0(adata.dataset$obs$ann,"_",adata.dataset$obs$sample)
  
  write_h5ad(adata.dataset,paste0("data/HLCA/datasets/",x,"/sc_adata.h5ad"))
}
)

remove(adata)
gc()

tf.split <- Sys.time()
tf.split - t0.split
```


