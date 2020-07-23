

```python
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
  # the file that will store the analysis results
```

    scanpy==1.4.4.post1 anndata==0.6.22.post1 umap==0.3.8 numpy==1.16.2 scipy==1.3.1 pandas==0.23.4 scikit-learn==0.20.3 statsmodels==0.10.1 python-igraph==0.7.1 louvain==0.6.1



```python
sc.settings.set_figure_params(dpi=100,dpi_save=300)
```


```python
%config InlineBackend.print_figure_kwargs={'facecolor' : "w"}
```


```python
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
```

# creat adata


```python
adata = sc.read('raw_data.csv',cache=False).T
```

    --> This might be very slow. Consider passing `cache=True`, which enables much faster reading from a cache file.



```python
anno = pd.read_csv('anno_data.csv',index_col = 0)
```


```python
adata.obs = anno
```


```python
filter_cells = pd.read_csv('filter_cells.csv',header = None)[0]
```


```python
adata.obs['spikeration'] = ''
adata.obs.loc[adata.obs_names.isin(filter_cells),'spikeration'] = 'T'
adata.obs.loc[~adata.obs_names.isin(filter_cells),'spikeration'] = 'F'
```


```python
adata = adata[adata.obs['Status'].isin(['Pre-target','Post-target'])]
```


```python
adata_raw = adata.copy()
```


```python
adata_raw.obs['n_counts'] = adata_raw.X.sum(axis=1)
adata_raw.obs['n_genes'] = (adata_raw.X > 0).sum(axis=1)
```


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

    filtered out 10588 genes that are detectedin less than 5 cells


    Trying to set attribute `.var` of view, making a copy.



```python
adata = adata[adata.obs['spikeration'] == 'F']
```


```python
sc.pp.filter_genes(adata, min_cells=5)
```

    filtered out 4 genes that are detectedin less than 5 cells


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




    False    14407
    True      1319
    Name: highly_variable, dtype: int64




```python
sc.pl.highly_variable_genes(adata)
```


![png](output_26_0.png)



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
sc.pp.regress_out(adata, 'Batch')
```

    regressing out ['n_genes', 'n_counts']


    ... storing 'Batch' as categorical
    ... storing 'Barcode' as categorical
    ... storing 'location' as categorical
    ... storing 'type' as categorical
    ... storing 'sex' as categorical
    ... storing 'sample' as categorical
    ... storing 'Clinical' as categorical
    ... storing 'Status' as categorical
    ... storing 'spikeration' as categorical


        finished (0:00:34)
    regressing out Batch
        finished (0:00:38)



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
sc.pl.pca_variance_ratio(adata,log = False)
```


![png](output_36_0.png)


## Computing the neighborhood graph 


```python
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=10,random_state = 0)
```

    computing neighbors
        using 'X_pca' with n_pcs = 10
        finished: added to `.uns['neighbors']`
        'distances', distances for each pair of neighbors
        'connectivities', weighted adjacency matrix (0:00:06)


## Embedding the neighborhood graph


```python
sc.tl.tsne(adata,n_pcs=10,random_state=5,perplexity = 40)
sc.pl.tsne(adata,size = 30)
```

    computing tSNE
        using 'X_pca' with n_pcs = 10
        using the 'MulticoreTSNE' package by Ulyanov (2017)
        finished: added
        'X_tsne', tSNE coordinates (adata.obsm) (0:00:07)



![png](output_40_1.png)


## Clustering the neighborhood graph 


```python
sc.tl.louvain(adata,resolution=0.9)
```

    running Louvain clustering
        using the "louvain" package of Traag (2017)
        finished: found 8 clusters and added
        'louvain', the cluster labels (adata.obs, categorical) (0:00:00)



```python
sc.pl.tsne(adata, color=['louvain','Status','sample'],ncols = 2,size = 40)
```


![png](output_43_0.png)



```python
adata.obs['cluster'] = ""
adata.obs.loc[adata.obs['louvain']=="0",'cluster'] = 'CD16+ Mo'
adata.obs.loc[adata.obs['louvain']=='2','cluster'] = 'pDC'
adata.obs.loc[adata.obs['louvain']=='1','cluster'] = 'CD14+ Mo'
adata.obs.loc[adata.obs['louvain']=="3",'cluster'] = 'cDC'
adata.obs.loc[adata.obs['louvain']=="6",'cluster'] = 'Plasma cell'
adata.obs.loc[adata.obs['louvain']=="5",'cluster'] = 'CD34+ Progenitor'
adata.obs.loc[adata.obs['louvain']=="4",'cluster'] = 'iMo'
adata.obs.loc[adata.obs['louvain']=="7",'cluster'] = 'T cell'
```


```python
sc.pl.tsne(adata, color='cluster',legend_loc = 'on data',legend_fontsize= 8)
```

    ... storing 'cluster' as categorical



![png](output_45_1.png)


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

    filtered out 13370 genes that are detectedin less than 3 cells



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


```python
adata_cdc
```




    AnnData object with n_obs × n_vars = 178 × 12948 
        obs: 'Batch', 'Barcode', 'location', 'type', 'sex', 'sample', 'Clinical', 'Status', 'cfBRAFV600E', 'spikeration', 'n_counts', 'n_genes'
        var: 'n_cells'



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




    False    11155
    True      1793
    Name: highly_variable, dtype: int64




```python
sc.pl.highly_variable_genes(adata_cdc)
```


![png](output_57_0.png)



```python
hvg = pd.Series(adata_cdc.var_names[adata_cdc.var['highly_variable']].tolist())
```


```python
set(hvg).intersection(set(sexgenes))
```




    {'RPS4Y1', 'TSIX', 'TXLNGY', 'ZFY'}




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




    False    11159
    True      1789
    Name: highly_variable, dtype: int64



# Regress_out & scale  (only genes & counts) 

## Regress_out


```python
sc.pp.regress_out(adata_cdc, ['n_genes','n_counts'])
sc.pp.regress_out(adata_cdc,'Batch')
```

    regressing out ['n_genes', 'n_counts']


    ... storing 'Batch' as categorical
    ... storing 'Barcode' as categorical
    ... storing 'location' as categorical
    ... storing 'type' as categorical
    ... storing 'sex' as categorical
    ... storing 'sample' as categorical
    ... storing 'Clinical' as categorical
    ... storing 'Status' as categorical
    ... storing 'spikeration' as categorical


        finished (0:00:21)
    regressing out Batch
        finished (0:00:23)



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
sc.pl.pca_variance_ratio(adata_cdc,log = False,n_pcs=50)
```


![png](output_68_0.png)



```python
sc.pl.pca_loadings(adata_cdc,components=[12,13,14,15])
```


![png](output_69_0.png)


## Computing the neighborhood graph 


```python
sc.pp.neighbors(adata_cdc, n_neighbors=8, n_pcs=12,random_state = 0)
```

    computing neighbors
        using 'X_pca' with n_pcs = 12
        finished: added to `.uns['neighbors']`
        'distances', distances for each pair of neighbors
        'connectivities', weighted adjacency matrix (0:00:00)



```python
sc.tl.umap(adata_cdc,min_dist=0.1,random_state=3)
sc.pl.umap(adata_cdc,size = 80)
```

    computing UMAP
        finished: added
        'X_umap', UMAP coordinates (adata.obsm) (0:00:02)



![png](output_72_1.png)


## Clustering the neighborhood graph 


```python
sc.tl.louvain(adata_cdc,resolution=0.4)
```

    running Louvain clustering
        using the "louvain" package of Traag (2017)
        finished: found 2 clusters and added
        'louvain', the cluster labels (adata.obs, categorical) (0:00:00)



```python
sc.pl.umap(adata_cdc, color=['louvain'],ncols = 2,size = 100)
```


![png](output_75_0.png)



```python
sc.pl.umap(adata_cdc, color=['CD1C','HLA-DQA1','HLA-DQB1'],ncols = 3,size = 100,cmap = 'PuRd')
```


![png](output_76_0.png)



```python
sc.pl.umap(adata_cdc, color=['CLEC9A','CADM1','THBD'],ncols = 4,size = 100,cmap = 'PuRd')
```


![png](output_77_0.png)



```python
adata_cdc.obs['cluster_cdc'] = adata_cdc.obs['louvain'].copy()
```


```python
new_cluster_names = ['DC2/DC3','DC1']
adata_cdc.rename_categories('cluster_cdc', new_cluster_names)
```


```python
orde = ['DC1','DC2/DC3']
adata_cdc.obs['cluster_cdc'] = pd.Categorical(adata_cdc.obs['cluster_cdc'],categories=orde,ordered=True)
```


```python
adata_cdc.uns['cluster_mdc_colors'] = ['#c02739',"#ffd082"]
sc.pl.umap(adata_cdc,color = ['cluster_cdc'],legend_loc = 'on data',size = 80)
```


![png](output_81_0.png)


## cluster_detail


```python
adata.obs['cluster_detail'] = list(adata.obs['cluster'].copy())
```


```python
adata.obs.loc[adata.obs_names.isin(adata_cdc.obs_names[adata_cdc.obs['cluster_cdc'] == 'DC1']),'cluster_detail'] = 'DC1'
adata.obs.loc[adata.obs_names.isin(adata_cdc.obs_names[adata_cdc.obs['cluster_cdc'] == 'DC2/DC3']),'cluster_detail'] = 'DC2/DC3'
```


```python
orde = ['CD14+ Mo','CD16+ Mo','iMo','DC1','DC2/DC3','pDC','T cell','Plasma cell','CD34+ Progenitor']
```


```python
adata.obs['cluster_detail'] = pd.Categorical(adata.obs['cluster_detail'],ordered=True,
                                          categories=orde)
```


```python
adata.uns['cluster_detail_colors'] = ['#e6739f',"#9818d6","#51eaea",'#c02739',"#ffd082","#2078b5","#4d4c7d",'#afc7e8',"#c1a57b"]
```


```python
sc.pl.tsne(adata,color = ['cluster_detail'],legend_fontsize = 8,title = '', size = 30,ncols =1,save = 'Cluster_detail.pdf')
```

    WARNING: saving figure to file figures/tsneCluster_detail.pdf



![png](output_88_1.png)


## degs_CD14+ Mo Pre_target vs. Post_target


```python
adata_14sub = adata[adata.obs['cluster'].isin(["CD14+ Mo"])]
adata_14sub
```




    View of AnnData object with n_obs × n_vars = 267 × 15726 
        obs: 'Batch', 'Barcode', 'location', 'type', 'sex', 'sample', 'Clinical', 'Status', 'cfBRAFV600E', 'spikeration', 'n_counts', 'n_genes', 'louvain', 'cluster', 'cluster_detail'
        var: 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
        uns: 'pca', 'neighbors', 'louvain', 'louvain_colors', 'Status_colors', 'sample_colors', 'cluster_colors', 'cluster_detail_colors'
        obsm: 'X_pca', 'X_tsne'
        varm: 'PCs'




```python
sc.tl.rank_genes_groups(adata_14sub, 'Status',method = 'wilcoxon',n_genes = 15726)
sc.tl.filter_rank_genes_groups(adata_14sub,
                               min_in_group_fraction=0.1,
                               max_out_group_fraction=1,
                               min_fold_change=0)
```

    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
    Filtering genes using: min_in_group_fraction: 0.1 min_fold_change: 0, max_out_group_fraction: 1


## degs_CD16+ Mo Pre_target vs. Post_target


```python
adata_16sub = adata[adata.obs['cluster'].isin(["CD16+ Mo"])]
adata_16sub
```




    View of AnnData object with n_obs × n_vars = 274 × 15726 
        obs: 'Batch', 'Barcode', 'location', 'type', 'sex', 'sample', 'Clinical', 'Status', 'cfBRAFV600E', 'spikeration', 'n_counts', 'n_genes', 'louvain', 'cluster', 'cluster_detail'
        var: 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
        uns: 'pca', 'neighbors', 'louvain', 'louvain_colors', 'Status_colors', 'sample_colors', 'cluster_colors', 'cluster_detail_colors'
        obsm: 'X_pca', 'X_tsne'
        varm: 'PCs'




```python
sc.tl.rank_genes_groups(adata_16sub, 'Status',method = 'wilcoxon',n_genes = 15726)
sc.tl.filter_rank_genes_groups(adata_16sub,
                               min_in_group_fraction=0.1,
                               max_out_group_fraction=1,
                               min_fold_change=0.0)
```

    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
    Filtering genes using: min_in_group_fraction: 0.1 min_fold_change: 0.0, max_out_group_fraction: 1


## degs_DC2/DC3 Pre_target vs. Post_target


```python
adata_dc2dc3 = adata[adata.obs['cluster_detail'].isin(['DC2/DC3'])]
adata_dc2dc3
```




    View of AnnData object with n_obs × n_vars = 166 × 15726 
        obs: 'Batch', 'Barcode', 'location', 'type', 'sex', 'sample', 'Clinical', 'Status', 'cfBRAFV600E', 'spikeration', 'n_counts', 'n_genes', 'louvain', 'cluster', 'cluster_detail'
        var: 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
        uns: 'pca', 'neighbors', 'louvain', 'louvain_colors', 'Status_colors', 'sample_colors', 'cluster_colors', 'cluster_detail_colors'
        obsm: 'X_pca', 'X_tsne'
        varm: 'PCs'




```python
sc.tl.rank_genes_groups(adata_dc2dc3, 'Status',method = 'wilcoxon',n_genes = 15726)
sc.tl.filter_rank_genes_groups(adata_dc2dc3,
                               min_in_group_fraction=0.1,
                               max_out_group_fraction=1,
                               min_fold_change=0)
```

    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
    Filtering genes using: min_in_group_fraction: 0.1 min_fold_change: 0, max_out_group_fraction: 1


## degs_iMo Pre_target vs. Post_target


```python
adata_iMo = adata[adata.obs['cluster'].isin(["iMo"])]
adata_iMo
```




    View of AnnData object with n_obs × n_vars = 131 × 15726 
        obs: 'Batch', 'Barcode', 'location', 'type', 'sex', 'sample', 'Clinical', 'Status', 'cfBRAFV600E', 'spikeration', 'n_counts', 'n_genes', 'louvain', 'cluster', 'cluster_detail'
        var: 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
        uns: 'pca', 'neighbors', 'louvain', 'louvain_colors', 'Status_colors', 'sample_colors', 'cluster_colors', 'cluster_detail_colors'
        obsm: 'X_pca', 'X_tsne'
        varm: 'PCs'




```python
sc.tl.rank_genes_groups(adata_iMo, 'Status',method = 'wilcoxon',n_genes = 15726)
sc.tl.filter_rank_genes_groups(adata_iMo,
                               min_in_group_fraction=0.1,
                               max_out_group_fraction=1,
                               min_fold_change=0.0)
```

    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
    Filtering genes using: min_in_group_fraction: 0.1 min_fold_change: 0.0, max_out_group_fraction: 1


## degs_pDC Pre_target vs. Post_target


```python
adata_pdc = adata[adata.obs['cluster'].isin(["pDC"])]
adata_pdc
```




    View of AnnData object with n_obs × n_vars = 252 × 15726 
        obs: 'Batch', 'Barcode', 'location', 'type', 'sex', 'sample', 'Clinical', 'Status', 'cfBRAFV600E', 'spikeration', 'n_counts', 'n_genes', 'louvain', 'cluster', 'cluster_detail'
        var: 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
        uns: 'pca', 'neighbors', 'louvain', 'louvain_colors', 'Status_colors', 'sample_colors', 'cluster_colors', 'cluster_detail_colors'
        obsm: 'X_pca', 'X_tsne'
        varm: 'PCs'




```python
sc.tl.rank_genes_groups(adata_pdc, 'Status',method = 'wilcoxon',n_genes = 15726)
sc.tl.filter_rank_genes_groups(adata_pdc,
                               min_in_group_fraction=0.1,
                               max_out_group_fraction=1,
                               min_fold_change=0)
```

    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
    Filtering genes using: min_in_group_fraction: 0.1 min_fold_change: 0, max_out_group_fraction: 1

