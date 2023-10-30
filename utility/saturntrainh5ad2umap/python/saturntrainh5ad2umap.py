#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import pandas as pd
import scanpy as sc
import seaborn as sns

x=sc.read(infile)
sc.tl.pca(x, random_state=seed)
sc.pp.neighbors(x, random_state=seed, n_neighbors=neighbor, n_pcs=npc)
sc.tl.leiden(x, resolution=resolution, random_state=seed)
sc.tl.umap(x, random_state=seed)

sc.set_figure_params(dpi_save=500, figsize=(5, 5))
for splitby in list(label)+['leiden']:
	ncolor=len(x.obs[splitby].value_counts())
	if ncolor<100:
		sc.pl.umap(x, color=splitby, frameon=False, show=False, save=f'{bname}_umap_{splitby}_wolabel.png')
		sc.pl.umap(
			x, color=splitby, frameon=False, show=False, save=f'{bname}_umap_{splitby}_ondata.png',
			legend_loc='on data', legend_fontsize='xx-small', legend_fontweight='normal'
			)
	else:
		palette=sns.husl_palette(ncolor)
		sc.pl.umap(x, color=splitby, palette=palette, frameon=False, show=False, save=f'{bname}_umap_{splitby}_wolabel.png')
		sc.pl.umap(x, color=splitby, palette=palette, frameon=False, show=False, save=f'{bname}_umap_{splitby}_ondata.png') ## duplicate

sc.write(filename=f'{bname}.h5ad', adata=x)
x.obs['barcode']=x.obs.index
x.obs.to_csv(f'{bname}_obs.txt.gz', sep='\t', index=False)
x.var['symbol']=x.var.index
x.var.to_csv(f'{bname}_var.txt.gz', sep='\t', index=False)
