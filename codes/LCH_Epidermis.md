

```python
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
```

    scanpy==1.4.4.post1 anndata==0.6.22.post1 umap==0.3.8 numpy==1.16.2 scipy==1.3.1 pandas==0.23.4 scikit-learn==0.20.3 statsmodels==0.10.1 python-igraph==0.7.1 louvain==0.6.1



```python
import matplotlib.pyplot as plt
```


```python
sc.settings.set_figure_params(dpi=100,dpi_save=300)
```


```python
%config InlineBackend.print_figure_kwargs={'facecolor' : "w"}
```


```python
plt.rcParams['pdf.fonttype'] = 42
```


```python
def degToTable(degs):
    degs_names =pd.DataFrame(degs['names']).melt(var_name='group', value_name='gene')
    degs_pvals_adj = pd.DataFrame(degs['pvals_adj']).melt(var_name='group', value_name='pvals_adj')
    degs_pvals = pd.DataFrame(degs['pvals']).melt(var_name='group', value_name='pvals')
    degs_logF = pd.DataFrame(degs['logfoldchanges']).melt(var_name='group', value_name='logfoldchanges')
    degs_merge = pd.concat([degs_names,degs_pvals,degs_pvals_adj,degs_logF],axis = 1)
    degs_merge = degs_merge.T.drop_duplicates().T
    return degs_merge
```

# creat adata


```python
adata_skin = sc.read('raw_data.csv',cache=False).T
```

    --> This might be very slow. Consider passing `cache=True`, which enables much faster reading from a cache file.



```python
anno_skin = pd.read_csv('anno_data.csv',index_col = 0)
```


```python
adata_skin.obs = anno_skin
```


```python
filter_cells = pd.read_csv('filter_cells.csv',header = None)[0]
```


```python
adata_skin.obs['spikeration'] = ''
adata_skin.obs.loc[adata_skin.obs_names.isin(filter_cells),'spikeration'] = 'T'
adata_skin.obs.loc[~adata_skin.obs_names.isin(filter_cells),'spikeration'] = 'F'
```


```python
set(adata_skin.obs['location'])
```




    {'Normal skin', 'PB_LCH', 'PB_control', 'Skin_LCH'}




```python
adata_skin = adata_skin[adata_skin.obs['location'].isin(['Skin_LCH','Normal skin'])]
```


```python
adata_skin_raw = adata_skin.copy()
```


```python
adata_skin_raw.obs['n_counts'] = adata_skin_raw.X.sum(axis=1)
adata_skin_raw.obs['n_genes'] = (adata_skin_raw.X > 0).sum(axis=1)
```


```python
sc.pl.scatter(adata_skin_raw, x='n_counts', y='n_genes')
```

    ... storing 'Batch' as categorical
    ... storing 'Barcode' as categorical
    ... storing 'location' as categorical
    ... storing 'type' as categorical
    ... storing 'sex' as categorical
    ... storing 'sample' as categorical
    ... storing 'Clinical' as categorical
    ... storing 'Status' as categorical
    ... storing 'spikeration' as categorical



![png](output_16_1.png)


# Basic filter


```python
adata_skin = adata_skin_raw[adata_skin_raw.obs_names[adata_skin_raw.obs['spikeration'] == 'F']].copy()
```


```python
adata_skin = adata_skin[(adata_skin.obs['n_genes'] > 1000),:]
adata_skin = adata_skin[(adata_skin.obs['n_counts'] > 1000)  & (adata_skin.obs['n_counts'] < 1000000), :]
```


```python
sc.pp.filter_genes(adata_skin,min_cells=5)
```

    filtered out 11991 genes that are detectedin less than 5 cells


    Trying to set attribute `.var` of view, making a copy.


# normalize (library-size correct)


```python
adata_skin.obs['n_counts'].describe()
```




    count       403.000000
    mean     105943.976562
    std      120364.921875
    min        2956.000000
    25%       19141.000000
    50%       53181.000000
    75%      152817.000000
    max      673898.000000
    Name: n_counts, dtype: float64




```python
sc.pp.normalize_per_cell(adata_skin, counts_per_cell_after=1e4)
```

    normalizing by total count per cell
        finished (0:00:00): normalized adata.X and added    'n_counts', counts per cell before normalization (adata.obs)



```python
sc.pp.log1p(adata_skin)
```


```python
adata_skin.raw = adata_skin
```

# Identify highly-variable genes.


```python
sc.pp.highly_variable_genes(adata_skin, min_mean=0.15,max_mean=10, min_disp=1)
```

    extracting highly variable genes
        finished (0:00:00)
    --> added
        'highly_variable', boolean vector (adata.var)
        'means', float vector (adata.var)
        'dispersions', float vector (adata.var)
        'dispersions_norm', float vector (adata.var)



```python
adata_skin.var['highly_variable'].value_counts()
```




    False    13223
    True      1104
    Name: highly_variable, dtype: int64




```python
sc.pl.highly_variable_genes(adata_skin)
```


![png](output_29_0.png)



```python
hvg = pd.Series(adata_skin.var_names[adata_skin.var['highly_variable']].tolist())
```

## Regress_out


```python
sc.pp.regress_out(adata_skin, ['n_genes','n_counts'])
```

    regressing out ['n_genes', 'n_counts']
        finished (0:00:24)



```python
sc.pp.regress_out(adata_skin,['Batch'])
```

    regressing out ['Batch']
        finished (0:00:24)



```python
sc.pp.scale(adata_skin, max_value=10)
```

## Principal component analysis


```python
sc.tl.pca(adata_skin, svd_solver='arpack',n_comps = 50)
```

    computing PCA with n_comps = 50
    computing PCA on highly variable genes
        finished (0:00:00)



```python
sc.pl.pca_variance_ratio(adata_skin,log =False,n_pcs=30)
```


![png](output_37_0.png)


## Computing the neighborhood graph 


```python
sc.pp.neighbors(adata_skin, n_neighbors=15, n_pcs=4,random_state = 0)
```

    computing neighbors
        using 'X_pca' with n_pcs = 4
        finished: added to `.uns['neighbors']`
        'distances', distances for each pair of neighbors
        'connectivities', weighted adjacency matrix (0:00:05)


## Embedding the neighborhood graph


```python
sc.tl.umap(adata_skin,min_dist= 0.5,random_state=6)
```

    computing UMAP
        finished: added
        'X_umap', UMAP coordinates (adata.obsm) (0:00:02)



```python
sc.tl.louvain(adata_skin,resolution=0.05)
```

    running Louvain clustering
        using the "louvain" package of Traag (2017)
        finished: found 2 clusters and added
        'louvain', the cluster labels (adata.obs, categorical) (0:00:00)



```python
sc.pl.umap(adata_skin,color=['louvain','location','sample','Batch','FCGR3B','CD3D','CD1C','CD1A','CD207','PTPRC','HLA-DRB1'],ncols = 4,size = 100,frameon = False)
```


![png](output_43_0.png)



```python
adata_skin_1 = adata_skin.copy()
```


```python
adata_skin_1.obs['cluster_raw'] = adata_skin_1.obs['louvain'].copy()
```


```python
adata_skin_1.obs['cluster_raw'] = list(adata_skin_1.obs['cluster_raw'])
adata_skin_1.obs.loc[adata_skin_1.obs['cluster_raw'] == '0','cluster_raw'] = 'Non-hemaptopoietic cell'
adata_skin_1.obs.loc[adata_skin_1.obs['cluster_raw'] == '1','cluster_raw'] = 'Hematopoetic cell'
```


```python
adata_skin = adata_skin_1[adata_skin_1.obs['louvain'] == '0']
```

# adata_skin


```python
adata_skin = adata_skin_raw[adata_skin.obs_names]
```


```python
sc.pp.filter_genes(adata_skin,min_cells=5)
```

    filtered out 12451 genes that are detectedin less than 5 cells


    Trying to set attribute `.var` of view, making a copy.


# normalize (library-size correct)


```python
adata_skin.obs['n_counts'].describe()
```




    count       353.000000
    mean     110024.679688
    std      123890.914062
    min        2956.000000
    25%       19284.000000
    50%       58290.000000
    75%      160456.000000
    max      673898.000000
    Name: n_counts, dtype: float64




```python
sc.pp.normalize_per_cell(adata_skin, counts_per_cell_after=1e4)
```

    normalizing by total count per cell
        finished (0:00:00): normalized adata.X and added    'n_counts', counts per cell before normalization (adata.obs)



```python
sc.pp.log1p(adata_skin)
```


```python
adata_skin.raw = adata_skin
```

# Identify highly-variable genes.


```python
sc.pp.highly_variable_genes(adata_skin, min_mean=0.15,max_mean=10, min_disp=0.8)
```

    extracting highly variable genes
        finished (0:00:00)
    --> added
        'highly_variable', boolean vector (adata.var)
        'means', float vector (adata.var)
        'dispersions', float vector (adata.var)
        'dispersions_norm', float vector (adata.var)



```python
adata_skin.var['highly_variable'].value_counts()
```




    False    12450
    True      1417
    Name: highly_variable, dtype: int64




```python
sc.pl.highly_variable_genes(adata_skin)
```


![png](output_59_0.png)



```python
hvg = pd.Series(adata_skin.var_names[adata_skin.var['highly_variable']].tolist())
```


```python
sexgenes = pd.read_csv('~/database/sexgene.csv')['x'].tolist()
```


```python
set(hvg).intersection(set(sexgenes))
```




    {'PRKY', 'UTY'}




```python
adata_skin.var['highly_variable'][set(hvg).intersection(set(sexgenes))] = False
```

    /home/Public/BioSoft/anaconda2/envs/python36/lib/python3.6/site-packages/ipykernel_launcher.py:1: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      """Entry point for launching an IPython kernel.


## Regress_out


```python
sc.pp.regress_out(adata_skin, ['n_genes','n_counts'])
```

    regressing out ['n_genes', 'n_counts']
        finished (0:00:25)



```python
sc.pp.regress_out(adata_skin,['Batch'])
```

    regressing out ['Batch']
        finished (0:00:25)



```python
sc.pp.scale(adata_skin, max_value=10)
```

## Principal component analysis


```python
sc.tl.pca(adata_skin, svd_solver='arpack',n_comps = 50)
```

    computing PCA with n_comps = 50
    computing PCA on highly variable genes
        finished (0:00:00)



```python
sc.pl.pca_variance_ratio(adata_skin,log =False,n_pcs=30)
```


![png](output_70_0.png)


## Computing the neighborhood graph 


```python
sc.pp.neighbors(adata_skin, n_neighbors=12, n_pcs=9,random_state = 10)
```

    computing neighbors
        using 'X_pca' with n_pcs = 9
        finished: added to `.uns['neighbors']`
        'distances', distances for each pair of neighbors
        'connectivities', weighted adjacency matrix (0:00:00)


## Embedding the neighborhood graph


```python
sc.tl.umap(adata_skin,min_dist= 0.4,random_state=38)
```

    computing UMAP
        finished: added
        'X_umap', UMAP coordinates (adata.obsm) (0:00:00)



```python
sc.tl.louvain(adata_skin,resolution=0.6)
sc.pl.umap(adata_skin,color=['louvain','sample','location','Batch'],ncols = 2,size = 80,frameon = False)
```

    running Louvain clustering
        using the "louvain" package of Traag (2017)
        finished: found 8 clusters and added
        'louvain', the cluster labels (adata.obs, categorical) (0:00:00)



![png](output_75_1.png)



```python
sc.tl.rank_genes_groups(adata_skin,groupby='louvain',method = 'wilcoxon',n_genes= 5000)
```

    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)



```python
sc.tl.filter_rank_genes_groups(adata_skin,min_in_group_fraction=0.1,max_out_group_fraction=0.8,min_fold_change=0.5)
```

    Filtering genes using: min_in_group_fraction: 0.1 min_fold_change: 0.5, max_out_group_fraction: 0.8



```python
degs_filter = adata_skin.uns['rank_genes_groups_filtered']
degs_filter = degToTable(degs_filter)
degs_filter = degs_filter.dropna()
degs_filter = degs_filter[degs_filter['pvals_adj'] < 0.05]
```


```python
adata_skin.obs['cluster_skin'] = adata_skin.obs['louvain'].copy() 
```


```python
new_cluster_names = ['LC-like1','LC-like2','LC','Prolif. LC','T cell','Neu','Mac','pDC']
adata_skin.rename_categories('cluster_skin', new_cluster_names)
```

# names clusters


```python
adata_skin.obs['cluster_skin'] = pd.Categorical(adata_skin.obs['cluster_skin'],categories=['LC','LC-like1','LC-like2','Prolif. LC','T cell','Mac','Neu','pDC'])
```


```python
adata_skin.uns['cluster_skin_colors'] = ['#14868c','#ff7315','#fe346e','#be7b98',"#81d4f9",'#933e15',"#85a392",'#152b93']
```


```python
sc.pl.umap(adata_skin,color=['cluster_skin'],ncols = 1,size = 100)
```


![png](output_84_0.png)



```python
sc.tl.rank_genes_groups(adata_skin,groupby='cluster_skin',method = 'wilcoxon',n_genes= 5000)
```

    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)



```python
sc.tl.filter_rank_genes_groups(adata_skin,min_in_group_fraction=0.1,max_out_group_fraction=0.8,min_fold_change=0.5)
```

    Filtering genes using: min_in_group_fraction: 0.1 min_fold_change: 0.5, max_out_group_fraction: 0.8

