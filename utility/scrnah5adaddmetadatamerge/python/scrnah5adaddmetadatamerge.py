#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

if __name__ == '__main__':
	import scanpy as sc
	import sys
	if f.endswith('.h5ad'):
		x=sc.read(f)
	else:
		print(f'Error: .h5ad file is expected. {f}')
		sys.exit(-1)

	import pandas as pd
	metadata=pd.read_table(metadata, header=0)
	print(metadata)
	print(metadata.shape)
	print(metadata.dtypes)

	tmp=x.obs
	tmp['_tmpindex_']=tmp.index
	if tmp[key].dtype!=metadata[key].dtype: # bug fix for incompatible dtye of the key, such as int64 vs character
		tmp[key]=tmp[key].astype(str)
		metadata[key]=metadata[key].astype(str)
	tmp=tmp.merge(metadata, how='inner', on=key)
	tmp.set_index('_tmpindex_', inplace=True)
	tmp.rename_axis(index=None, inplace=True)

	x=x[tmp.index].copy() # to ensure order between x.X and x.obs
	x.obs=tmp
	if len(character)>0:
		for h in character:
			if h in x.obs:
				x.obs[h]=x.obs[h].astype(str)
	print(x.obs)
	print(x.obs.shape)
	print(x.obs.dtypes)

	sc.write(filename=f'{outdir}/{bname}.h5ad', adata=x)
	x.obs['barcode']=x.obs.index
	x.obs.to_csv(f'{outdir}/{bname}_obs.txt.gz', sep='\t', index=False)
	x.var['symbol']=x.var.index
	x.var.to_csv(f'{outdir}/{bname}_var.txt.gz', sep='\t', index=False)
