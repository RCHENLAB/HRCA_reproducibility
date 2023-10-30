#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad

if infile.endswith('.h5ad'):
	x=sc.read(infile)
else:
	print(f'Error: format is not supported for {infile}')
	sys.exit(-1)
print('==> x')
print(x)

if filterlabel:
	x=x[x.obs.groupby(filterlabel).filter(lambda x: len(x)>=mincell).index].copy()
	print('==> x, after filtering')
	print(x)

def colSums(obs):
	return np.sum(x[obs.index, ].X, axis=0)

xx=x.obs.groupby(key).apply(colSums)
print('==> xx')
print(xx)
print(xx.index)
print(xx.values)

tmp=ad.AnnData(
	X=np.asarray(np.stack(xx.values)),
	obs=pd.DataFrame(
		xx.index.tolist(),
		index=map(lambda x: separator.join(x), xx.index.values) if len(key)>1 else xx.index.values,
		columns=key,
		),
	var=pd.DataFrame(index=x.var.index),
	)
tmp.obs.insert(loc=0, column=idheader, value=tmp.obs.index)
tmp.var.insert(loc=0, column='symbol', value=tmp.var.index)
print('==> tmp')
print(tmp)
print('==> tmp.X')
print(tmp.X)
print('==> tmp.obs')
print(tmp.obs)
print('==> tmp.var')
print(tmp.var)

sc.write(filename=f'{outdir}/{bname}.h5ad', adata=tmp)
tmp.obs.to_csv(f'{outdir}/{bname}_obs.txt.gz', sep='\t', index=False)
tmp.var.to_csv(f'{outdir}/{bname}_var.txt.gz', sep='\t', index=False)
