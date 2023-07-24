--- 
title: "Metacell Tutorial"
author: "Mariia Bilous, Léonard Hérault, Aurélie Gabriel, David Gfeller"
date: "2023-07-24"
site: bookdown::bookdown_site
documentclass: book
bibliography:
- packages.bib
description: |
  This is a tutorial about metacell construction and analysis.
link-citations: yes
github-repo: GfellerLab/Metacell_tutorial
---

# About

The structure of this tutorial.





## Book structure

Book consists of several Chapters (i.e., first-level headings). Each chapter is in separate .Rmd file in the root folder, with a name in XY_text.Rmd format, with `XY` being numbers.

Each Chapter consists of sections and sup-sections (i.e., second-level and lower heading), files for which are located in `./sub_pages`.  Sub-pages and chapters may also call *functional_chunks*, which are located in `./functional_chunks` and represent parts of code that can be repetitively run (e.g., `load_anndata`, `save_mc_anndata` etc). When the book is rendered, the included *sub_pages* and *functional_chunks* are basically inserted in the Chapter as inline code. The only challenge is the relative path of the files and resulting outputs, such as plots. To resolve this issue, currently, I manually set up the project folder as a knitting root directory in each sub-file (i.e., sub_pages and functional_chunks) as `knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())` . Also, in my RStudio settings, I have the following setting `Tools -> Global Options... -> R Markdown -> Evaluate chunk in directory -> Project`.

**Note:** each chapters runs in a new R session and they do not share the environment, thus, we need to provide global knit options for each chapter, otherwise they are lost. I do it with a `source('./R/config.R')` in the beginning of each chapter. 


## Installation and requirements

R requirements

```r
install.packages('rprojroot') # to reset work directory to the Project root
install.packages('bookdown') # to render book
```

To run **MC2** and **SEACells** in RStudio, we need 

```r
install.packages('reticulate') # to run Python
```

Then, we need to setup virtual environment


```bash
pip install virtualenv
cd <Path_to_Metacell_tutorial>
virtualenv my_env
source my_env/bin/activate

# Installing SEACells, pip install installs old version, that does not work for me, thus install from git
git clone https://github.com/dpeerlab/SEACells.git
cd SEACells
python setup.py install

cd ..

pip install -r SEACells_requirements.txt # here some packages have wrong/non-existing vision, so I manually changed their versions 
pip install ipywidgets
pip install jupyter

pip install metacells

# in project dir
echo 'RETICULATE_PYTHON=my_env/bin/python' > '.Renviron' 

# restart RStudio and open 'Metacell_tutorial.Rproj' 

```

## Render book 

The function to render book is `bookdown::render_book()`, this will take some time, as it will execute all the chunks in the book, there is an option to cache some chunks, but we have to make sure that cached chunks do not share variables with non-cached chunks (it will raise an error anyway). 

`bookdown::preview_chapter()` renders a chapter.

## Get data


To get 3k PBMCs, use scanpy datasets af follows 

```python
import scanpy as sc 
import os

ad = sc.datasets.pbmc3k()
adata_proc = sc.datasets.pbmc3k_processed()

ad       = ad[adata_proc.obs_names].copy()
ad.obs   = adata_proc.obs.copy()
ad.uns   = adata_proc.uns.copy()
ad.obsm  = adata_proc.obsm.copy()
ad.obsp  = adata_proc.obsp.copy()

raw_ad = sc.AnnData(ad.X.copy())
raw_ad.obs_names, raw_ad.var_names = ad.obs_names, ad.var_names
ad.raw = raw_ad

sc.pl.embedding(ad, 'X_umap', color='louvain')
#> /mnt/c/Aurelie/postdoc_UNIL/Metacell_tutorial/my_env/lib/python3.8/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
#>   cax = scatter(
```

<img src="index_files/figure-html/unnamed-chunk-5-1.png" width="672" />
<!-- We pre-process the data using standard steps from scanpy .... -->
<!-- ```{r supercell-standard-preproc, child='./functional_chunks/scanpy_standard_preproc.Rmd'} -->
<!-- ``` -->



```python
directory = os.path.join("data", "3k_pbmc")

if not os.path.exists(directory):
    os.makedirs(directory)
    
ad.write_h5ad(os.path.join("data", "3k_pbmc", "singlecell_anndata_filtered.h5ad"))
```



```r
library(reticulate)
library(Seurat)
#> Attaching SeuratObject
raw_counts <- Matrix::t(as(py$ad$raw$X, "CsparseMatrix"))
colnames(raw_counts) <- rownames(py$ad$obs)
rownames(raw_counts) <- rownames(py$ad$var)

# norm_counts <- Matrix::t(as(py$ad$X, "CsparseMatrix"))
# colnames(norm_counts) <- rownames(py$ad$obs)
# rownames(norm_counts) <- rownames(py$ad$var)

pbmc <- CreateSeuratObject(counts = raw_counts, meta.data = py$ad$obs)
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes
#> ('-')
# pbmc@assays$RNA@data <- norm_counts
saveRDS(pbmc, file = paste0("data/3k_pbmc/singlecell_seurat_filtered.rds"))
```
