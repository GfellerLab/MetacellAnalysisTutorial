import pandas as pd
import numpy as np
import palantir
import warnings

# adapted from https://github.com/dpeerlab/SEACells/blob/main/SEACells/evaluate.py
def celltype_frac(x, col_name):
 val_counts = x[col_name].value_counts()
 return val_counts.values[0] / val_counts.values.sum()
 
# adapted from https://github.com/dpeerlab/SEACells/blob/main/SEACells/evaluate.py
def purity(ad, annotation_label, MC_label = 'membership'):
	celltype_fraction = ad.obs.groupby(MC_label).apply(lambda x: celltype_frac(x, annotation_label))
	celltype = ad.obs.groupby(MC_label).apply(lambda x: x[annotation_label].value_counts().index[0])
	
	return pd.concat([celltype, celltype_fraction], axis=1).rename(columns={0: annotation_label, 1: f'{annotation_label}_purity'})
	

# adapted from https://github.com/dpeerlab/SEACells/blob/main/SEACells/evaluate.py

def compactness(ad, low_dim_embedding = 'X_pca', MC_label = 'SEACell', DO_DC = True, name = 'compactness', n_comp = None):
    """
    Compute compactness of each metacell. Compactness is defined is the average variance of diffusion components 
    across cells that constitute a metcell

    :param ad: (Anndata) Anndata object
    :param low_dim_embedding: (str) `ad.obsm` field for constructing diffusion components
    :param SEACell_label: (str) `ad.obs` field for computing diffusion component variances

    :return: `pd.DataFrame` with a dataframe of compactness per metacell

    """
    
    components = pd.DataFrame(ad.obsm[low_dim_embedding]).set_index(ad.obs_names)
    
    if n_comp is None:
        n_comp = components.shape[1]
    else:
        if n_comp > components.shape[1]:
            warning(f'Number of components to use is larges than number of existing components, n_comp set to max')
            n_comp = components.shape[1]
    
    components = components.iloc[:,:n_comp]
    
    if DO_DC:
        dm_res = palantir.utils.run_diffusion_maps(components,  n_components=n_comp)
        emb = palantir.utils.determine_multiscale_space(dm_res, n_eigs=n_comp)
    else:
        emb = components

    res = pd.DataFrame(emb.join(ad.obs[MC_label]).groupby(MC_label).var().mean(1)).rename(columns={0:name})
    res['low_dim_embedding'] = low_dim_embedding
    res['n_comp'] = n_comp
    res[MC_label] = res.index.to_list()
    
    if DO_DC:
        res['low_dim_embedding'] = 'DC'
    
    return res
  
  
# adapted from https://github.com/dpeerlab/SEACells/blob/main/SEACells/evaluate.py
def separation(
    ad,
    low_dim_embedding='X_pca',
    DO_DC = True,
    n_comp = None,
    name = 'separation',
    nth_nbr=1,
    cluster=None,
    MC_label='SEACell'):
    """
    Compute separation of each metacell. Separation is defined is the distance to the nearest neighboring metacell

    :param ad: (Anndata) Anndata object
    :param low_dim_embedding: (str) `ad.obsm` field for constructing diffusion components
    :param nth_nbr: (int) Which neighbor to use for computing separation
    :param SEACell_label: (str) `ad.obs` field for computing diffusion component variances

    :return: `pd.DataFrame` with a separation of compactness per metacell

    """
    components = pd.DataFrame(ad.obsm[low_dim_embedding]).set_index(ad.obs_names)
    
    if n_comp is None:
        n_comp = components.shape[1]
    else:
        if n_comp > components.shape[1]:
            warning(f'Number of components to use is larges than number of existing components, n_comp set to max')
            n_comp = components.shape[1]
    
    components = components.iloc[:,:n_comp]
    
    if DO_DC:
        dm_res = palantir.utils.run_diffusion_maps(components, n_components=max(10, n_comp))
        dc = palantir.utils.determine_multiscale_space(dm_res, n_eigs=n_comp)
    else:
        dc = components
        

    # Compute DC per metacell
    metacells_dcs = dc.join(ad.obs[MC_label], how='inner').groupby(MC_label).mean()

    from sklearn.neighbors import NearestNeighbors
    neigh = NearestNeighbors(n_neighbors=nth_nbr)
    nbrs = neigh.fit(metacells_dcs)
    dists, nbrs = nbrs.kneighbors()
    dists = pd.DataFrame(dists).set_index(metacells_dcs.index)
    dists.columns += 1

    nbr_cells = np.array(metacells_dcs.index)[nbrs]

    metacells_nbrs = pd.DataFrame(nbr_cells)
    metacells_nbrs.index = metacells_dcs.index
    metacells_nbrs.columns += 1

    res = pd.DataFrame(dists[nth_nbr]).rename(columns={1:name})

    res['low_dim_embedding'] = low_dim_embedding
    res['n_comp'] = n_comp
    res[MC_label] = res.index.to_list()
    
    if DO_DC:
        res['low_dim_embedding'] = 'DC'
        
    return res

# adapted from https://github.com/tanaylab/metacells/blob/master/metacells/tools/quality.py 
import scipy

def mc_gene_var(
	ad, 
	MC_label,
	):
		
		"""
		Get mc gene variance
		"""
		if scipy.sparse.issparse(ad.X):
			X = pd.DataFrame(ad.X.A)
		else:
			X = pd.DataFrame(ad.X)
			
		X = X.groupby(ad.obs[MC_label].to_list()).var()
		X = X.set_axis(ad.var_names, axis=1, inplace=False)
		return X

def mc_gene_mean(
	ad, 
	MC_label
	):
		
	"""
	Get mc gene variance
	"""
	if scipy.sparse.issparse(ad.X):
	 	X = pd.DataFrame(ad.X.A)
	else:
		X = pd.DataFrame(ad.X)
	
	X = X.groupby(ad.obs[MC_label].to_list()).mean()
	X = X.set_axis(ad.var_names, axis=1, inplace=False)
	return X


def mc_inner_normalized_var(
		ad,
		MC_label
):
	"""
	Gene normalized variance within metacells
	"""
	v = mc_gene_var(ad, MC_label)
	m = mc_gene_mean(ad, MC_label)
	zeros_mask = m == 0
		
	res = np.reciprocal(m, where = zeros_mask)
	res[zeros_mask] = 0
	res *= v
		
	res[zeros_mask] = np.nan
	return res
		
	
