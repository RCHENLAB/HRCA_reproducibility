#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import seaborn as sns
import scanpy as sc

x=sc.read(f)

if 'X_umap' not in x.obsm:
	print('Error: X_umap is missing. See scanpy.tl.umap()')
	sys.exit(-1)

sc.set_figure_params(dpi_save=500, figsize=(width, height))

if len(label)==0:
	label=x.obs.columns

for splitby in label:
	if splitby not in x.obs_keys():
		print(f'Error: {splitby} is not a metadata.')
		continue
		# sys.exit(-1)

	ncolor=len(x.obs[splitby].value_counts())
	if ncolor<100:
		sc.pl.umap(x, color=splitby, frameon=False, show=False, title=title, save=f"{bname}_umap_{splitby}_wolabel.{format}")
		sc.pl.umap(x, color=splitby, frameon=False, show=False, title=title, save=f"{bname}_umap_{splitby}_ondata.{format}",
			legend_loc='on data',
			legend_fontsize='xx-small',
			legend_fontweight='normal',
			)
		sc.pl.umap(x, color=splitby, frameon=False, show=False, title=title, save=f"{bname}_umap_{splitby}_fontline.{format}",
			legend_loc='on data',
			legend_fontsize='xx-small',
			legend_fontweight='normal',
			legend_fontoutline=1,
			)
	else:
		palette=sns.husl_palette(ncolor)
		sc.pl.umap(x, color=splitby, palette=palette, frameon=False, show=False, title=title, save=f"{bname}_umap_{splitby}_wolabel.{format}")
		sc.pl.umap(x, color=splitby, palette=palette, frameon=False, show=False, title=title, save=f"{bname}_umap_{splitby}_ondata.{format}") # duplicate
