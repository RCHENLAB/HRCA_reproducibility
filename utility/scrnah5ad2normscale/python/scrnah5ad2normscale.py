#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
if f.endswith('.h5ad'):
	x=sc.read(f)
elif f.endswith('.h5'):
	x=sc.read_10x_h5(f)
else:
	print(f'Error: format is not supported for {f}')
	import sys
	sys.exit(-1)

if useraw:
	x=x.raw.to_adata()
elif raw:
	x.raw=x

sc.pp.normalize_total(x, target_sum=targetsum)
sc.pp.log1p(x)

if scale:
	sc.pp.scale(x)

sc.write(filename=f'{outdir}/{bname}.h5ad', adata=x)
x.obs['barcode']=x.obs.index
x.obs.to_csv(f'{outdir}/{bname}_obs.txt.gz', sep='\t', index=False)
x.var['symbol']=x.var.index
x.var.to_csv(f'{outdir}/{bname}_var.txt.gz', sep='\t', index=False)
