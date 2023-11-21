# Command line tool for MC construction
Add before that a QC section: and say that in a benchmark setting, some parameters need to be the same: k, nb pca, gamma, filtering?, latent space.
In the next section we present command line tools that enables users to test the 3 tools to build metacells.

Here is a command line tool to construct metacells using either tool (MC2, SuperCell or SEACells) from a provided dataset.


```bash
setwd("MetacellToolkit_path")
```


```bash
proj_name="3k_pbmc"
MC_tool="SuperCell"
# input raw adata output adata
Rscript cli/${MC_tool}CL.R -i data/${proj_name}/singlecell_anndata_filtered.h5ad -o data/${proj_name}/${MC_tool}/ -n 50 -f 2000 -k 30 -g 50 -s adata

# input raw adata output seurat
Rscript cli/${MC_tool}CL.R -i data/${proj_name}/singlecell_anndata_filtered.h5ad -o data/${proj_name}/${MC_tool}/ -n 50 -f 2000 -k 30 -g 50 -s seurat
#> $ARGS
#> character(0)
#> 
#> $input
#> [1] "data/3k_pbmc/singlecell_anndata_filtered.h5ad"
#> 
#> $outdir
#> [1] "data/3k_pbmc/SuperCell/"
#> 
#> $n.pc
#> [1] 50
#> 
#> $nFeatures
#> [1] 2000
#> 
#> $k.knn
#> [1] 30
#> 
#> $gamma
#> [1] 50
#> 
#> $output
#> [1] "adata"
#> 
#> $nPCs
#> [1] 50
#> 
#> $isNorm
#> [1] FALSE
#> 
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#>           used  (Mb) gc trigger  (Mb) max used  (Mb)
#> Ncells 2694980 144.0    4872347 260.3  4872347 260.3
#> Vcells 7938289  60.6   18438552 140.7 14602449 111.5
#> [1] "Normalize data..."
#> [1] "Identify Metacells..."
#> [1] "Identify Metacells sequentially..."
#> [1] "Treat SeuratProject cells"
#> [1] "Identify 53 metacells using SuperCell..."
#> Warning: The following arguments are not used: row.names
#> [1] "Assign metadata to metacells and compute purities..."
#> $ARGS
#> character(0)
#> 
#> $input
#> [1] "data/3k_pbmc/singlecell_anndata_filtered.h5ad"
#> 
#> $outdir
#> [1] "data/3k_pbmc/SuperCell/"
#> 
#> $n.pc
#> [1] 50
#> 
#> $nFeatures
#> [1] 2000
#> 
#> $k.knn
#> [1] 30
#> 
#> $gamma
#> [1] 50
#> 
#> $output
#> [1] "seurat"
#> 
#> $nPCs
#> [1] 50
#> 
#> $isNorm
#> [1] FALSE
#> 
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#>           used  (Mb) gc trigger  (Mb) max used  (Mb)
#> Ncells 2694980 144.0    4872518 260.3  4872518 260.3
#> Vcells 7938289  60.6   18438552 140.7 14602450 111.5
#> [1] "Normalize data..."
#> [1] "Identify Metacells..."
#> [1] "Identify Metacells sequentially..."
#> [1] "Treat SeuratProject cells"
#> [1] "Identify 53 metacells using SuperCell..."
#> Warning: The following arguments are not used: row.names
#> [1] "Assign metadata to metacells and compute purities..."
```



```bash
proj_name="3k_pbmc"
MC_tool="SEACells"
# input raw adata output adata
Rscript cli/${MC_tool}CL.R -i data/${proj_name}/singlecell_anndata_filtered.h5ad -o data/${proj_name}/${MC_tool}/ -n 50 -f 2000 -k 30 -g 50 -s adata

# input raw adata output seurat
Rscript cli/${MC_tool}CL.R -i data/${proj_name}/singlecell_anndata_filtered.h5ad -o data/${proj_name}/${MC_tool}/ -n 50 -f 2000 -k 30 -g 50 -s seurat
```
