#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
import sys
if f.endswith('.h5ad'):
	tmp=sc.read(f)
else:
	print(f'Error: .h5ad file is expected. {f}')
	sys.exit(-1)

import numpy as np
vs=np.unique(tmp.obs[splitby])
for v in vs:
	x=tmp[tmp.obs[splitby]==v].copy()
	sc.write(filename=f'{outdir}/{bname}_{v}.h5ad', adata=x)
	x.obs['barcode']=x.obs.index
	x.obs.to_csv(f'{outdir}/{bname}_{v}_obs.txt.gz', sep='\t', index=False)
	x.var['symbol']=x.var.index
	x.var.to_csv(f'{outdir}/{bname}_{v}_var.txt.gz', sep='\t', index=False)
