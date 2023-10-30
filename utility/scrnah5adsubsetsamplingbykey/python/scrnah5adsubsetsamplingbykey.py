#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
if f.endswith('.h5ad'):
	x=sc.read(f)
else:
	print(f'Error: format is not supported for {f}')
	import sys
	sys.exit(-1)

print('\n\nBefore sampling:')
print(x.X.shape)
print(x.obs[key].value_counts())

tmp=x.obs[key].value_counts()>nsample # categories with sufficent number of cells for sampling
s1=x.obs[x.obs[key].isin(tmp[tmp].index)]
s2=x.obs[~x.obs[key].isin(tmp[tmp].index)]
if s1.shape[0]>0:
	barcode=s1.groupby(key, observed=True).sample(n=nsample, random_state=seed).index.to_list()+s2.index.to_list()
else:
	barcode=s2.index.to_list()

x=x[barcode].copy()
print('\n\nAfter sampling:')
print(x.X.shape)
print(x.obs[key].value_counts())

sc.write(filename=f'{outdir}/{bname}.h5ad', adata=x)
x.obs['barcode']=x.obs.index
x.obs.to_csv(f'{outdir}/{bname}_obs.txt.gz', sep='\t', index=False)
x.var['symbol']=x.var.index
x.var.to_csv(f'{outdir}/{bname}_var.txt.gz', sep='\t', index=False)
