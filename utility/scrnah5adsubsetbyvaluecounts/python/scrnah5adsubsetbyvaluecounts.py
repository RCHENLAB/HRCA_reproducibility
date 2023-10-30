#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

if __name__ == '__main__':
	import scanpy as sc
	if f.endswith('.h5ad'):
		x=sc.read(f)
	else:
		print("Error: please input .h5ad file.")
		import sys
		sys.exit(-1)
	import pandas as pd
	xx=pd.value_counts(x.obs[label])>=ncell
	x=x[x.obs[label].isin(xx[xx].index)].copy()
	sc.write(filename=f'{bname}.h5ad', adata=x)
	x.obs['barcode']=x.obs.index
	x.obs.to_csv(f'{bname}_obs.txt.gz', sep='\t', index=False)
	x.var['symbol']=x.var.index
	x.var.to_csv(f'{bname}_var.txt.gz', sep='\t', index=False)
