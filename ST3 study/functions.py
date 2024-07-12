# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:41:49 2019

@author: dmorgan
"""
import os
import numpy as np
import pandas as pd
import anndata as an
import scanpy as sc
import igraph
import statsmodels.api as sm
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import sklearn.preprocessing
import pickle
from statsmodels.formula.api import ols
import statsmodels

def subset(origdata, cells):
    df = origdata.raw.X[origdata.obs_names.isin(cells),:]
    adata = an.AnnData(df)
    adata.obs_names = origdata.obs_names[origdata.obs_names.isin(cells)]
    adata.var_names = origdata.raw.var_names
    adata.obs = origdata.obs.loc[adata.obs_names,:]
    return adata

def process(adata, model = 'linear', npcs = 10):
    adata.raw = adata
    sc.pp.filter_genes_dispersion(adata, flavor = 'seurat')
    if (model == 'linear'):
        sc.pp.regress_out(adata, keys = 'n_counts')
        sc.pp.scale(adata)
        print('linear scaling')
    else:
        adata = poissonscale(adata, 'n_counts')
        print('poisson scaling')
    sc.tl.pca(adata, svd_solver = 'arpack')
    sc.pp.neighbors(adata, n_pcs = npcs, n_neighbors = 30)
    sc.tl.umap(adata, min_dist = .3)
    sc.tl.leiden(adata, resolution = None)
    sc.pl.umap(adata,color ='leiden')
    return adata



def getgenes(adata, cluster, field = 'leiden', thresh = .25):
    data1 = adata.raw.X[adata.obs[field].isin(cluster),:]
    data2 = adata.raw.X[~adata.obs[field].isin(cluster),:]
    genes = adata.raw.var_names
    g1sums = (data1>0).sum(axis = 0)
    g2sums = (data2>0).sum(axis = 0)
    genes_use = genes[g1sums + g2sums > (data2.shape[1] + data1.shape[1])*.03]
    mean1 = np.log(np.mean(np.expm1(data1), axis = 0) + 1)
    mean2 = np.log(np.mean(np.expm1(data2), axis = 0) + 1)
    logfc = mean1 - mean2
    pval = scipy.stats.ttest_ind(data1, data2, axis = 0, equal_var = False)
    pct_1 = ((data1 > 0).sum(axis = 0)) / data1.shape[0]
    pct_2 = ((data2 > 0).sum(axis = 0)) / data2.shape[0]
    df = pd.DataFrame({'LogFC' : logfc,'pct1' : pct_1, 'pct2' : pct_2,  'pval' : pval.pvalue})
    df.index = genes
    df = df.sort_values(by = 'pval', ascending = True)
    df = df[(df.pct1 > .1) | (df.pct2 > .1)]
    df = df[abs(df.LogFC) > thresh]
    df['pval_adj'] = statsmodels.stats.multitest.multipletests(df.pval.values, method = 'bonferroni')[1]
    return df.sort_values(by = 'LogFC', ascending = False)

def addmeta(adata, meta):
    adata.obs['tumor'] = pd.Categorical(meta.loc[adata.obs_names, 'Tumor'])
    adata.obs['size'] = pd.Categorical(meta.loc[adata.obs_names, 'Size'])
    adata.obs['source'] = pd.Categorical(meta.loc[adata.obs_names, 'Source'])
    return adata

def umap(adata, color, size = 10):
    coords = pd.DataFrame(adata.obsm.X_umap, columns = ['x', 'y'])
    coords['color'] = adata.obs.loc[:, color].values
    coords = coords.sample(frac = 1)
    sns.lmplot(x = 'x', y = 'y', data = coords, hue = 'color', fit_reg = False, palette = "Set1", scatter_kws={"s": size})

def diffgenes(adata, c1, c2, field = 'leiden', threshold = .25):
    data1 = adata.raw.X[adata.obs[field].isin(c1),:]
    data2 = adata.raw.X[adata.obs[field].isin(c2),:]
    genes = adata.raw.var_names
    g1sums = (data1>0).sum(axis = 0)
    g2sums = (data2>0).sum(axis = 0)
    genes_use = genes[g1sums + g2sums > (data2.shape[1] + data1.shape[1])*.03]
    mean1 = np.log(np.mean(np.expm1(data1), axis = 0) + 1)
    mean2 = np.log(np.mean(np.expm1(data2), axis = 0) + 1)
    logfc = mean1 - mean2
    pval = scipy.stats.ttest_ind(data1, data2, equal_var = False)
    pct_1 = ((data1 > 0).sum(axis = 0)) / data1.shape[0]
    pct_2 = ((data2 > 0).sum(axis = 0)) / data2.shape[0]
    df = pd.DataFrame({'LogFC' : logfc,'pct1' : pct_1, 'pct2' : pct_2,  'pval' : pval.pvalue})
    df.index = genes
    df = df.sort_values(by = 'pval', ascending = True)
    df = df[(df.pct1 > .1) | (df.pct2 > .1)]
    df = df[abs(df.LogFC) > threshold]
    return df.sort_values(by = 'LogFC', ascending = False)

def subprocess(adata, cells, model = 'linear'):
    adata = subset(adata, cells)
    adata = process(adata, model)
    return adata

def removecluster(adata, cluster, field = 'leiden', model = 'linear'):
    obs = adata.obs_names[~adata.obs[field].isin(cluster)]
    return subprocess(adata, obs,  model)

def whichcells(adata, cluster, field = 'leiden'):
     return adata.obs_names[adata.obs[field].isin(cluster)]   
    

    
def diffusion(adata, root, field = 'leiden', nd = 2):
    adata.uns['iroot'] = np.flatnonzero(adata.obs[field] == root)[0]
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata, n_dcs = nd)
    sc.pl.diffmap(adata, color=['dpt_pseudotime'])
    return adata

def dplot(adata, color = 'dpt_pseudotime', dim1 = 1, dim2 = 2):
    coords = pd.DataFrame(adata.obsm.X_diffmap[:,[(dim1), (dim2)]] , columns = ['x', 'y'])
    coords['color'] = adata.X[:,adata.var_names == color]
    coords = coords.sample(frac = 1)
    plt.scatter(x = coords['x'], y = coords['y'], c = coords['color'])
    

def barplot(tab, index1 = 'size', index2 = 'leiden'):
    tab = tab.groupby([index1, index2]).size().unstack()
    tab = tab.fillna(0)
    tab = tab.div(tab.sum(axis = 1), axis = 0)
    if (index1 == 'size'):
        tab = tab.reindex(['S', 'M', "L", 'X'])
    elif (index1 == 'tumor'):
        inds = ['S1' ,"S3", 'S4', 'S6', 
               'M1', 'M2', 'M3', 'M4', 
               'L1', 'L2', 'L3', 'L4',
               'X1', 'X2', 'X3', 'X4']
        inds = [x for x in inds if (x in set(tab.index))]
        tab = tab.reindex(inds)
    tab.plot.bar(stacked = True).legend(bbox_to_anchor =(1,1))
    return tab

def pseudogene(adata, gene):
    from patsy import dmatrix
    import statsmodels.api as sm
    time = adata.obs.dpt_pseudotime.values
    gene = adata.raw.X[:,adata.raw.var_names == gene]
    df = pd.DataFrame({'x' : time, 'y' : gene[:,0]})
    t_x = dmatrix("bs(df.x, df = 3)", {'df.x' : df.x}, return_type = 'dataframe')
    fit = sm.GLM(df.y, t_x).fit()
    pred = fit.predict(dmatrix('bs(df.x, df = 3)', {'df.x' : df.x}, return_type = 'dataframe'))
    plt.plot(time, gene, 'o', time, pred, 'o')
    plt.show()
    
def scaleRow(row, regressOut):
    import statsmodels.api as sm
    poisson_model = sm.GLM(row, regressOut, family = sm.families.Poisson())
    results = poisson_model.fit()
    adjusted = results.resid_pearson - min(results.resid_pearson)
    scaled = sklearn.preprocessing.scale(adjusted, axis = 0)
    return scaled

def poissonscale(adata, regressOut):
    mat = adata.X
    regressOut = adata.obs[regressOut]
    for i in range(0, adata.X.shape[1]):
        row = mat[:,i]
        regressOut = sm.add_constant(regressOut)
        mat[:,i] = scaleRow(row, regressOut)
    adata.X = mat
    return adata

def summaryplots(adata): 
    umap(adata, 'leiden')
    umap(adata, 'source')
    umap(adata, 'size')
    
def barsource(adata, index1 = 'size', index2 = 'leiden'):
    a = barplot(adata.obs[adata.obs.source == "Tumo"], index1, index2)
    a = barplot(adata.obs[adata.obs.source == "CD45"], index1, index2)
    a = barplot(adata.obs[adata.obs.source == "PBMC"], index1, index2)
    
def cellnova(tab,  index):
    tab = tab.groupby(['tumor', index]).size().unstack()
    tab = tab.fillna(0)
    tab = tab.div(tab.sum(axis = 1), axis = 0)
    names = tab.index.astype(str)
    sizes = [a[0] for a in names]
    
    for curr in tab.columns:
        print("\n" + curr)
        df = pd.DataFrame({'Size': sizes, 'Var': tab[curr]})
        anova = ols('Var ~ Size', data = df).fit()
        anova = sm.stats.anova_lm(anova, typ = 2)
        print(anova)
        
        
def seuratExport(adata, fname):
    df = pd.DataFrame(adata.raw.X.transpose())
    df.columns = adata.obs_names
    df.index = adata.raw.var_names
    df = df.reset_index()
    pyarrow.feather.write_dataframe(df, fname + '.feather') 
    adata.obs.to_csv(fname + '_meta.txt')
    
def importSeurat(filename):
    rawdata = pd.read_feather(filename + '_raw.feather')
    metadata = pd.read_feather(filename + '_meta.feather')
    metadata.index = metadata.cell
    rawdata.index = rawdata.gene.values
    rawdata = rawdata.drop('gene', axis = 1)
    adata = an.AnnData(rawdata.transpose())
    sc.pp.filter_genes(adata, min_cells = 3)
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    adata = process(adata, npcs = 10)
    adata.obs = adata.obs.join(metadata, rsuffix = 'r')
    return adata