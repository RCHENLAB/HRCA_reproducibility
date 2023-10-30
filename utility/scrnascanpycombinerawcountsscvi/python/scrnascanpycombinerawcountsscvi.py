#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
tmp=sc.read(scvi)

obs=tmp.obs
if len(obss)>0:
	setobs=set(obss)
	if invert:
		obss=[k for k in tmp.obs_keys() if k not in setobs]
	else:
		obss=[k for k in tmp.obs_keys() if k in setobs]
	obs=tmp.obs[obss]
print(obs)

x=sc.read(rawcount)
x=x[obs.index].copy()

import anndata as ad
x=ad.AnnData(X=x.X, obs=obs, var=x.var)
x.obsm['X_scVI']=tmp.obsm['X_scVI']
x.obsm['X_umap']=tmp.obsm['X_umap']

sc.write(filename=f'{outdir}/{bname}.h5ad', adata=x)
x.obs['barcode']=x.obs.index
x.obs.to_csv(f'{outdir}/{bname}_obs.txt.gz', sep='\t', index=False)
x.var['symbol']=x.var.index
x.var.to_csv(f'{outdir}/{bname}_var.txt.gz', sep='\t', index=False)
