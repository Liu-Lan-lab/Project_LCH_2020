

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
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
```


```python
sc.settings.set_figure_params(dpi=100,dpi_save=300)
```


```python
%config InlineBackend.print_figure_kwargs={'facecolor' : "w"}
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

# LCH_Epidermis & LCH_TreatmentNaive_PB


```python
adata_skin = sc.read('LCH_TreatmentNaive_PB.h5ad')
adata_pb = sc.read('LCH_Epidermis.h5ad')
```


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
adata.obs['n_counts'] = adata.X.sum(axis=1)
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
```


```python
adata = adata[adata.obs_names.isin(set(adata_skin.obs_names).union(adata_pb.obs_names))]
```


```python
sc.pp.filter_genes(adata,min_cells=5)
```

    filtered out 8001 genes that are detectedin less than 5 cells


    Trying to set attribute `.var` of view, making a copy.



```python
adata.obs['cluster_raw'] = ''

adata.obs.loc[adata.obs_names.isin(adata_skin.obs_names[adata_skin.obs['cluster_skin'] == 'LC-like1']),'cluster_raw'] = 'LC-like1'
adata.obs.loc[adata.obs_names.isin(adata_skin.obs_names[adata_skin.obs['cluster_skin'] == 'LC-like2']),'cluster_raw'] = 'LC-like2'
adata.obs.loc[adata.obs_names.isin(adata_skin.obs_names[adata_skin.obs['cluster_skin'] == 'LC']),'cluster_raw'] = 'LC'
adata.obs.loc[adata.obs_names.isin(adata_skin.obs_names[adata_skin.obs['cluster_skin'] == 'Prolif. LC']),'cluster_raw'] = 'Prolif. LC'
adata.obs.loc[adata.obs_names.isin(adata_skin.obs_names[adata_skin.obs['cluster_skin'] == 'Neu']),'cluster_raw'] = 'Neu_skin'
adata.obs.loc[adata.obs_names.isin(adata_skin.obs_names[adata_skin.obs['cluster_skin'] == 'Mac']),'cluster_raw'] = 'Mac_skin'
adata.obs.loc[adata.obs_names.isin(adata_skin.obs_names[adata_skin.obs['cluster_skin'] == 'T cell']),'cluster_raw'] = 'T cell_skin'
adata.obs.loc[adata.obs_names.isin(adata_skin.obs_names[adata_skin.obs['cluster_skin'] == 'pDC']),'cluster_raw'] = 'pDC_skin'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'B cell']),'cluster_raw'] = 'B cell'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'CD14+ Mo']),'cluster_raw'] = 'CD14+ Mo'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'CD16+ Mo']),'cluster_raw'] = 'CD16+ Mo'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'DC1']),'cluster_raw'] = 'DC1'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'DC2/DC3']),'cluster_raw'] = 'DC2/DC3'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'pDC']),'cluster_raw'] = 'pDC'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'CD34+ Progenitor']),'cluster_raw'] = 'CD34+ Progenitor'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'iMo']),'cluster_raw'] = 'iMo'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'Neu']),'cluster_raw'] = 'Neu'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'T cell']),'cluster_raw'] = 'T cell'
adata.obs.loc[adata.obs_names.isin(adata_pb.obs_names[adata_pb.obs['cluster_detail2'] == 'Plasma cell']),'cluster_raw'] = 'Plasma cell'
```


```python
adata_raw = adata.copy()
```

# normalize (library-size correct).


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
sc.pp.highly_variable_genes(adata, min_mean=0.5, max_mean=10, min_disp=1)
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




    False    17246
    True      1071
    Name: highly_variable, dtype: int64




```python
sc.pl.highly_variable_genes(adata)
```


![png](output_22_0.png)



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




```python
hvg_skin = pd.Series(list(adata_skin.var_names[adata_skin.var['highly_variable']]))
hvg_pb = pd.Series(list(adata_pb.var_names[adata_pb.var['highly_variable']]))
```


```python
hvg_merge = set(set(hvg_skin).union(set(hvg_pb)))
```


```python
adata.var['highly_variable'][hvg_merge] = True
```

    /home/Public/BioSoft/anaconda2/envs/python36/lib/python3.6/site-packages/ipykernel_launcher.py:1: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      """Entry point for launching an IPython kernel.



```python
adata.var['highly_variable'].value_counts()
```




    False    15852
    True      2465
    Name: highly_variable, dtype: int64



# Regress_out & scale


```python
sc.pp.regress_out(adata, ['n_genes','n_counts'])
sc.pp.regress_out(adata,'Batch')
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
    ... storing 'cluster_raw' as categorical


        finished (0:00:49)
    regressing out Batch
        finished (0:01:08)



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
plt.rcParams['figure.figsize'] = 5,5
```


```python
sc.pl.pca_variance_ratio(adata,log = False)
```


![png](output_36_0.png)


## Computing the neighborhood graph 


```python
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=9,random_state = 6)
```

    computing neighbors
        using 'X_pca' with n_pcs = 9
        finished: added to `.uns['neighbors']`
        'distances', distances for each pair of neighbors
        'connectivities', weighted adjacency matrix (0:00:06)



```python
sc.tl.umap(adata,min_dist=0.4,random_state=14,n_components=2)
```

    computing UMAP
        finished: added
        'X_umap', UMAP coordinates (adata.obsm) (0:00:08)



```python
sc.pl.umap(adata,color = 'cluster_raw',size = 20,components = '1,2')
```


![png](output_40_0.png)



```python
orde = ['CD14+ Mo','CD16+ Mo','iMo','DC1','DC2/DC3','pDC','Neu','T cell','B cell','Plasma cell','CD34+ Progenitor',
        'LC','LC-like1','LC-like2','MKI67+ LC-like','Neu_skin','T cell_skin','Mac_skin','pDC_skin']
```


```python
adata.obs['cluster_raw'] = pd.Categorical(adata.obs['cluster_raw'],ordered=True,categories=orde)
```


```python
adata.uns['cluster_raw_colors'] = ['#e6739f',"#9818d6","#51eaea",'#c02739',"#ffd082","#2078b5",'#937d14',"#4d4c7d","#00bd56",'#afc7e8',"#c1a57b",
                            '#14868c','#ff7315','#fe346e','#be7b98',"#152b93",'#81d4f9',"#933e15",'#b892f7']
```


```python
sc.pl.umap(adata, color=['cluster_raw'],size = 30,save  = 'cluster_merge.pdf')
```

    WARNING: saving figure to file figures/umapcluster_merge.pdf



![png](output_44_1.png)


# adata_sub : mDC,pDC，iMo，CD14,CD16,


```python
adata_sub = adata.copy()
```


```python
adata_sub = adata_sub[adata_sub.obs['cluster_raw'].isin(['CD14+ Mo','CD16+ Mo','DC1','DC2/DC3','pDC','LC','LC-like1','LC-like2','iMo','Prolif. LC','pDC_skin'])]
```


```python
plt.rcParams['figure.figsize'] = 5,5
```


```python
sc.pl.umap(adata_sub,color = 'cluster_raw',size = 40)
```


![png](output_49_0.png)


# adata_sub2: DC2/3, iMo，CD14,


```python
adata_sub2 = adata.copy()
```


```python
adata_sub2 = adata_sub2[adata_sub2.obs['cluster_raw'].isin(['CD14+ Mo','DC2/DC3','LC-like1','LC-like2','iMo','LC','DC1'])]
```


```python
adata_sub2.obs['cluster_raw'].value_counts()
```




    CD14+ Mo    474
    DC2/DC3     246
    iMo         206
    LC-like1    111
    LC-like2     76
    LC           61
    DC1          12
    Name: cluster_raw, dtype: int64



## Score


```python
markers_dc2dc3 = ["CD1C","CLEC10A", 'FCER1A', "HLA-DQA2", "CD1D", "ANXA1", "S100A9", "HLA-DQA1", "HLA-DPB1", "VCAN", "CD1B", "CD33", "CD200R1", "CD180","FCGR2B","KLF4", "SIGLEC1", "ITGAX","CD86", "IFI30", "CLEC4A"]
```


```python
sc.tl.score_genes(adata_sub2, gene_list=markers_dc2dc3,use_raw=True)
adata_sub2.obs['DC2DC3_score'] = adata_sub2.obs['score'].copy()
```

    computing score 'score'


    Trying to set attribute `.obs` of view, making a copy.


        finished: added
        'score', score of gene set (adata.obs) (0:00:00)



```python
markers_mono = ["CD14", "TREM1", "FCGR3A", "FCGR3B", "CD163", "TLR4", "CLEC7A", "TLR2", "ITGAM", "ITGB2", "CTSD", "CTSA", "BST1", "STAB1", "IRAK3", "NLRP3", "CD55", "CXCR1", "FCAR", "CD52", "C3AR", "SIGLEC10","S100A8"]
```


```python
sc.tl.score_genes(adata_sub2, gene_list=markers_mono,use_raw=True)
adata_sub2.obs['Mono_score'] = adata_sub2.obs['score'].copy()
```

    computing score 'score'
    WARNING: gene: C3AR is not in adata.var_names and will be ignored
        finished: added
        'score', score of gene set (adata.obs) (0:00:00)



```python
markers_dc1 = ["CLEC9A", "HLA-DPA1", "CADM1", "CAMK2D", "IDO1", "CLNK", "ZNF366", "NDRG2", "XCR1", "CD59", "SLAMF8", "CD141", "BTLA", "C1ORF54","HAVCR2"]
```


```python
sc.tl.score_genes(adata_sub2, gene_list=markers_dc1,use_raw=True)
adata_sub2.obs['dc1_score'] = adata_sub2.obs['score'].copy()
```

    computing score 'score'
    WARNING: gene: CD141 is not in adata.var_names and will be ignored
    WARNING: gene: C1ORF54 is not in adata.var_names and will be ignored
        finished: added
        'score', score of gene set (adata.obs) (0:00:00)



```python
obs_subdata = pd.DataFrame(adata_sub2.obs)
```
