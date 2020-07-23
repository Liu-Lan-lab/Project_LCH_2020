

```python
import numpy as np
import pandas as pd
import sys
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
  # the file that will store the analysis results
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
adata_raw = sc.read('raw_data.csv',cache=False).T
```

    --> This might be very slow. Consider passing `cache=True`, which enables much faster reading from a cache file.



```python
anno = pd.read_csv('anno_data.csv',index_col = 0)
```


```python
adata_raw.obs = anno
```


```python
sc.pp.filter_cells(adata_raw, min_genes=200)
sc.pp.filter_genes(adata_raw, min_cells=3)
```

    filtered out 41 cells that haveless than 200 genes expressed
    filtered out 5824 genes that are detectedin less than 3 cells



```python
filter_cells = pd.read_csv('filter_cells.csv',header = None)[0]
```


```python
adata_raw.obs['spikeration'] = ''
adata_raw.obs.loc[adata_raw.obs_names.isin(filter_cells),'spikeration'] = 'T'
adata_raw.obs.loc[~adata_raw.obs_names.isin(filter_cells),'spikeration'] = 'F'
```


```python
adata_raw = adata_raw[adata_raw.obs['location'].isin(['PB_LCH','PB_control'])]
```


```python
adata_raw = adata_raw[adata_raw.obs['Status'].isin(['LCH','Healthy children'])]
```


```python
adata_raw.obs['n_counts'] = adata_raw.X.sum(axis=1)
```

    Trying to set attribute `.obs` of view, making a copy.



```python
sc.pl.scatter(adata_raw, x='n_counts', y='n_genes',color = 'spikeration')
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



```python
cell_filter_1 = set((adata_raw.obs_names[adata_raw.obs['n_genes'] < 1000]) | (adata_raw.obs_names[adata_raw.obs['n_genes'] > 10000])).union(set(adata_raw.obs_names[adata_raw.obs['spikeration']=='T']))
len(cell_filter_1)
```




    103




```python
cell_filter_2 = set((adata_raw.obs_names[adata_raw.obs['n_counts'] < 1000]) | (adata_raw.obs_names[adata_raw.obs['n_counts'] > 1000000]))
len(cell_filter_2)
```




    12




```python
cell_filter = cell_filter_1.union(cell_filter_2)
len(cell_filter)
```




    105




```python
adata = adata_raw.copy()
```

# Basic filtering


```python
adata = adata[(adata.obs['n_genes'] > 1000) , :]
adata = adata[(adata.obs['n_counts'] > 1000) & (adata.obs['n_counts'] < 1000000), :]
```


```python
sc.pp.filter_genes(adata, min_cells=5)
```

    filtered out 2911 genes that are detectedin less than 5 cells


    Trying to set attribute `.var` of view, making a copy.



```python
adata = adata[adata.obs['spikeration'] == 'F']
```


```python
sc.pp.filter_genes(adata, min_cells=5)
```

    filtered out 11 genes that are detectedin less than 5 cells


    Trying to set attribute `.var` of view, making a copy.


# normalize (library-size correct)


```python
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e5)
```

    normalizing by total count per cell
        finished (0:00:00): normalized adata.X and added    'n_counts', counts per cell before normalization (adata.obs)



```python
sc.pp.log1p(adata)
```


```python
adata.raw = adata
```

# Identify highly-variable genes.


```python
sc.pp.highly_variable_genes(adata, min_mean=0.3, max_mean=10, min_disp=1)
```

    extracting highly variable genes
        finished (0:00:00)
    --> added
        'highly_variable', boolean vector (adata.var)
        'means', float vector (adata.var)
        'dispersions', float vector (adata.var)
        'dispersions_norm', float vector (adata.var)



```python
adata.var['highly_variable'].value_counts()
```




    False    16246
    True      1326
    Name: highly_variable, dtype: int64




```python
sc.pl.highly_variable_genes(adata)
```


![png](output_33_0.png)



```python
hvg = pd.Series(adata.var_names[adata.var['highly_variable']].tolist())
```


```python
sexgenes = pd.read_csv('~/database/sexgene.csv')['x'].tolist()
```


```python
set(hvg).intersection(set(sexgenes))
```




    set()



# Regress_out & scale  (only genes & counts) 

## Regress_out


```python
sc.pp.regress_out(adata, ['n_genes','n_counts'])
sc.pp.regress_out(adata,'Batch')
```

    regressing out ['n_genes', 'n_counts']
        finished (0:00:45)
    regressing out Batch
        finished (0:00:58)



```python
sc.pp.scale(adata, max_value=10)
```

## Principal component analysis


```python
sc.tl.pca(adata, svd_solver='arpack')
```

    computing PCA with n_comps = 50
    computing PCA on highly variable genes
        finished (0:00:00)



```python
sc.pl.pca_variance_ratio(adata,log = False,n_pcs= 40)
```


![png](output_43_0.png)


## Computing the neighborhood graph 


```python
sc.pp.neighbors(adata, n_neighbors=12,n_pcs=10,random_state = 0)
```

    computing neighbors
        using 'X_pca' with n_pcs = 10
        finished: added to `.uns['neighbors']`
        'distances', distances for each pair of neighbors
        'connectivities', weighted adjacency matrix (0:00:01)


## Embedding the neighborhood graph


```python
sc.tl.tsne(adata,n_pcs=17,random_state=6,perplexity = 50)
sc.pl.tsne(adata,size = 20)
```

    computing tSNE
        using 'X_pca' with n_pcs = 17
        using the 'MulticoreTSNE' package by Ulyanov (2017)
        finished: added
        'X_tsne', tSNE coordinates (adata.obsm) (0:00:14)


    ... storing 'cluster' as categorical



![png](output_47_2.png)


## Clustering the neighborhood graph 


```python
sc.tl.louvain(adata,resolution=0.9)
```

    running Louvain clustering
        using the "louvain" package of Traag (2017)
        finished: found 11 clusters and added
        'louvain', the cluster labels (adata.obs, categorical) (0:00:00)



```python
sc.pl.tsne(adata, color=['louvain','location','sample'],ncols = 2,size = 20)
```


![png](output_50_0.png)



```python
sc.pl.tsne(adata, color=['sex'],ncols = 1,size = 20)
```


![png](output_51_0.png)


### Names adata_clusters


```python
sc.pl.tsne(adata,color = ['louvain'],legend_loc = 'on data',size = 30)
```


![png](output_53_0.png)



```python
adata.obs['cluster'] = adata.obs['louvain'].copy()
```


```python
new_cluster_names = ['CD14+ Mo','CD16+ Mo','pDC','cDC','iMo','Plasma cell','16Mo','CD34+ Progenitor','Neu','B cell','T cell']
adata.rename_categories('cluster', new_cluster_names)
```


```python
adata.obs['cluster'] = list(adata.obs['cluster'])
adata.obs.loc[adata.obs['cluster'] == '16Mo','cluster'] = 'CD16+ Mo'
```

## cDC


```python
adata_cdc = adata[adata.obs['cluster'] == 'cDC']
```


```python
adata_cdc = adata_raw[adata_cdc.obs_names,:].copy()
```


```python
sc.pp.filter_genes(adata_cdc,min_cells=3)
```

    filtered out 6774 genes that are detectedin less than 3 cells



```python
adata_cdc
```




    AnnData object with n_obs × n_vars = 258 × 13720 
        obs: 'Batch', 'Barcode', 'location', 'type', 'sex', 'sample', 'Clinical', 'Status', 'cfBRAFV600E', 'n_genes', 'spikeration', 'n_counts'
        var: 'n_cells'
        uns: 'spikeration_colors'




```python
sc.pp.normalize_per_cell(adata_cdc, counts_per_cell_after=1e5)
```

    normalizing by total count per cell
        finished (0:00:00): normalized adata.X and added    'n_counts', counts per cell before normalization (adata.obs)



```python
sc.pp.log1p(adata_cdc)
```


```python
adata_cdc.raw = adata_cdc
```

# Identify highly-variable genes.


```python
sc.pp.highly_variable_genes(adata_cdc, min_mean=0.15,max_mean=10, min_disp=1)
```

    extracting highly variable genes
        finished (0:00:00)
    --> added
        'highly_variable', boolean vector (adata.var)
        'means', float vector (adata.var)
        'dispersions', float vector (adata.var)
        'dispersions_norm', float vector (adata.var)



```python
adata_cdc.var['highly_variable'].value_counts()
```




    False    11997
    True      1723
    Name: highly_variable, dtype: int64




```python
sc.pl.highly_variable_genes(adata_cdc)
```


![png](output_68_0.png)



```python
hvg = pd.Series(adata_cdc.var_names[adata_cdc.var['highly_variable']].tolist())
```


```python
set(hvg).intersection(set(sexgenes))
```




    {'LINC00278', 'USP9Y', 'XIST'}




```python
adata_cdc.var['highly_variable'][set(hvg).intersection(set(sexgenes))] = False
```

    /home/Public/BioSoft/anaconda2/envs/python36/lib/python3.6/site-packages/ipykernel_launcher.py:1: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      """Entry point for launching an IPython kernel.



```python
adata_cdc.var['highly_variable'].value_counts()
```




    False    12000
    True      1720
    Name: highly_variable, dtype: int64



# Regress_out & scale  (only genes & counts) 

## Regress_out


```python
sc.pp.regress_out(adata_cdc, ['n_genes','n_counts'])
sc.pp.regress_out(adata_cdc,'Batch')
```

    regressing out ['n_genes', 'n_counts']
        finished (0:00:23)
    regressing out Batch
        finished (0:00:28)



```python
sc.pp.scale(adata_cdc, max_value=10)
```

## Principal component analysis


```python
sc.tl.pca(adata_cdc, svd_solver='arpack')
```

    computing PCA with n_comps = 50
    computing PCA on highly variable genes
        finished (0:00:00)



```python
sc.pl.pca_variance_ratio(adata_cdc,log = False,n_pcs=30)
```


![png](output_79_0.png)


## Computing the neighborhood graph 


```python
sc.pp.neighbors(adata_cdc, n_neighbors=12, n_pcs=30,random_state = 0)
```

    computing neighbors
        using 'X_pca' with n_pcs = 30
        finished: added to `.uns['neighbors']`
        'distances', distances for each pair of neighbors
        'connectivities', weighted adjacency matrix (0:00:00)



```python
sc.tl.umap(adata_cdc,min_dist=0.2,random_state=6)
sc.pl.umap(adata_cdc,size = 80)
```

    computing UMAP
        finished: added
        'X_umap', UMAP coordinates (adata.obsm) (0:00:00)



![png](output_82_1.png)


## Clustering the neighborhood graph 


```python
sc.tl.louvain(adata_cdc,resolution=0.65)
```

    running Louvain clustering
        using the "louvain" package of Traag (2017)
        finished: found 2 clusters and added
        'louvain', the cluster labels (adata.obs, categorical) (0:00:00)



```python
sc.pl.umap(adata_cdc, color=['louvain'],ncols = 2,size = 60)
```


![png](output_85_0.png)



```python
sc.pl.umap(adata_cdc, color=['CD1C','HLA-DQA1','HLA-DQB1'],ncols = 3,size = 50,cmap = 'PuRd')
```


![png](output_86_0.png)



```python
sc.pl.umap(adata_cdc, color=['CLEC9A','XCR1','CADM1'],ncols = 3,size = 50,cmap = 'PuRd',save = 'DC1_feature.pdf')
```

    WARNING: saving figure to file figures/umapDC1_feature.pdf



![png](output_87_1.png)



```python
adata_cdc.obs['cluster_cdc'] = adata_cdc.obs['louvain'].copy()
```


```python
new_cluster_names = ['DC2/DC3','DC1']
adata_cdc.rename_categories('cluster_cdc', new_cluster_names)
```
