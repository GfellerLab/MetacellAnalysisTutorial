# Construction of metacells (for 'continuous' data) {#MC-continuous}





```
#> findfont: Font family ['Raleway'] not found. Falling back to DejaVu Sans.
#> findfont: Font family ['Lato'] not found. Falling back to DejaVu Sans.
```


## MC2 (Python)



[//]: # (Chunk to run MC-2 metacell construction for a continuous dataset)



**This is a template, not working code**




Here we construct metacells using [Metacell-2 (MC2)](https://github.com/tanaylab/metacells). The code is adapted from the [author's tutorial](https://metacells.readthedocs.io/en/latest/Metacells_Vignette.html).

#### Imports {-}

```python
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import metacells as mc

import sys
sys.path.append('./mc_QC/') 
import mc_QC
```


#### Parameters {-}

```python

MC_tool = "MC2"
gamma   = 50 # graining level

# Here we can modify dataset 
proj_name = "CD34"

annotation_label = "cell_type" # name of annotation field (obs)
```

#### Load data {-}



**Setup MC object**

```python
mc.ut.set_name(ad, proj_name)
```


#### Gene filtering according to MC2 {-}


```python
excluded_gene_names = [] # for example, ['IGHMBP2', 'IGLL1', 'IGLL5', 'IGLON5', 'NEAT1', 'TMSB10', 'TMSB4X']
excluded_gene_patterns = ['MT-.*']
mc.pl.analyze_clean_genes(ad,
                          excluded_gene_names=excluded_gene_names,
                          excluded_gene_patterns=excluded_gene_patterns,
                          random_seed=123456)

mc.pl.pick_clean_genes(ad)
```

#### Cell fintering {-}
Since our data is pre-processed and low-quality cells have been already filtered out, this step is not applicable to our data, but we keep this chunk so that you can apply cell filtering to a newly generated dataset.

The first round of cell cleaning usually implies filltering out cell with very low and very hight UMI content. The second round includes cell filtering based on mitochondrial and/or ribosomal content. We will skip both steps as our data have been pre-filtered and will use very lenient cutoffs (`properly_sampled_min_cell_total`, `properly_sampled_max_cell_total` and `properly_sampled_max_excluded_genes_fraction`) such that all the cells are kept for the metacell construction.


```python
### The first round (high/low UMIs)
properly_sampled_min_cell_total = {'cell_lines' : 5000, '3k_pbmc': 200}[proj_name] # setup for the dataset that will be used 
properly_sampled_max_cell_total = {'cell_lines' : 110000, '3k_pbmc': 10000}[proj_name] # setup for the dataset that will be used 
```





```python
## The second round (content of non-clean genes, e.g., mito-genes)
properly_sampled_max_excluded_genes_fraction = 0.25
```







### Running MC2 {-}

Metacell-2 uses its own feature selection approach (i.e., selection of genes used to build metacells). Additionally, we can explicitly specify which features to use by providing two arguments:  feature_gene_names - genes that have to be used  forbidden_gene_names - genes to exclude.
In contrast to the SuperCell and SEACells, Metacell-2 does not allow to explicitly obtain metacell data at a user-defined graining level. Instead, to vary graining level, we have to vary a `target_metacell_size` parameter, that is `160000` by default. Here we provide a chunk to calibrate this value to get a desired graining level. Please, increase or decrease scale if the obtained graining level `gamma_obtained` is lower or larger than the requested one (gamma).

**Estimate target_metacell_size (gamma)**

```python
print(f'The requested graining level is {gamma}, lets estimate the target_metacell_size that should result in such graining level.')

scale = 2 # increase or decrease if the obtained graining level (`gamma_obtained`) is significantly > or < then the requested one `gamma`

N_c = ad.shape[0]

# estimated mean UMI content in downsampled data
est_downsample_UMI = np.quantile(np.array(total_umis_of_cells), 0.05)

target_metacell_size = int(est_downsample_UMI * gamma * scale)
target_metacell_size
```

#### Aggregate metacells {.unnumbered #aggregate-mc2}

```python
mc.pl.divide_and_conquer_pipeline(
    ad,
    #feature_gene_names   = feature_gene_names, # comment this line to allow Metacell2 selecting features
    #forbidden_gene_names = forbidden_gene_names, # comment this line to allow Metacell2 selecting features
    target_metacell_size = target_metacell_size,
    random_seed = 123456)

## make anndata of metacells
mc_ad = mc.pl.collect_metacells(ad, name='metacells')
```

Here we estimate whether a deviation of the obtained gamma is acceptable, and if not, suggest to increase or decrease `scale` parameter to get better graining level.


If the obtained graining level is not acceptable and you updated `scale` parameter according to suggestion, do not forget to re-run chunk \@ref(chunk:mc2-divide-n-conquer) \@ref(aggregate-mc2)


#### Visualize metacells {-}

```python
mc.pl.compute_umap_by_features(mc_ad, max_top_feature_genes=1000,
                               min_dist=2.0, random_seed=123456)
                               
umap_x = mc.ut.get_o_numpy(mc_ad, 'umap_x')
umap_y = mc.ut.get_o_numpy(mc_ad, 'umap_y')
```


```python
plt.figure()
sns.scatterplot(x=umap_x, y=umap_y)
plt.show()
```


```python
# make a membership -- index of metacell single cell belongs to 
ad.obs['membership'] = [int(i)+1 if i >= 0 else np.nan for i in ad.obs.metacell] 

## Save single-cell metadata (i.e., `raw.obs` dataframe) in the metacell adata object
mc_ad.uns = ad.uns.copy()
mc_ad.uns['sc.obs'] = ad.obs.copy()

# save the requested gamma
mc_ad.uns['gamma'] = gamma
```

### Compute latent space for metacell QC {-}




### Metacell QC {-}









## SuperCell (R)




[//]: # (Code to run mc construction with SuperCell for a continuous dataset)

Under construction...



## SEACells (Python)



[//]: # (Code to run mc construction with SEACells for a continuous dataset)

Under construction...