#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
import scvi
import anndata
import seaborn as sns
import scipy

x=sc.read(f)
x.var_names_make_unique()
print(vars(x))

# bug fix filtering genes with zero count
gene_subset, _=sc.pp.filter_genes(x, min_counts=1, inplace=False)
x=x[:, gene_subset].copy()

## Delete because not used
## if csr:
## 	x.X=scipy.sparse.csr_matrix(x.X)

x.layers['hvgcounts']=x.X.copy()
sc.pp.normalize_total(x)
sc.pp.log1p(x)
x.raw=x

if len(x.var.index)>ntop:
	# seurat_v3 uses raw counts, while others use normalized counts
	if flavor=='seurat_v3':
		sc.pp.highly_variable_genes(
			x,
			flavor=flavor,
			n_top_genes=ntop,
			subset=True,
			layer='hvgcounts',
			batch_key=batchkey,
		)
	else:
		sc.pp.highly_variable_genes(
			x,
			flavor=flavor,
			n_top_genes=ntop,
			subset=True,
			batch_key=batchkey,
		)

scvi.model.SCVI.setup_anndata(x, layer='hvgcounts', batch_key=batchkey)
vae=scvi.model.SCVI(x, n_layers=nlayer, n_latent=nlatent)
vae.train(max_epochs=epoch)
vae.save(f'{bname}_model')
x.obsm['X_scVI']=vae.get_latent_representation()
if normcounts:
	denoised=anndata.AnnData(X=vae.get_normalized_expression(), obs=x.obs)
	sc.write(filename=f'{bname}_denoised.h5ad', adata=denoised)

sc.pp.neighbors(x, use_rep='X_scVI', random_state=seed)
sc.tl.leiden(x, resolution=0.618, random_state=seed)
sc.tl.umap(x, random_state=seed)

sc.set_figure_params(dpi_save=500, figsize=(5, 5))
for splitby in [batchkey, 'leiden']:
	ncolor=len(x.obs[splitby].value_counts())
	if ncolor<100:
		sc.pl.umap(x, color=splitby, frameon=False, show=False, save=f'{bname}_umap_{splitby}_wolabel.png')
		sc.pl.umap(x, color=splitby, frameon=False, show=False, save=f'{bname}_umap_{splitby}_ondata.png'
			, legend_loc='on data', legend_fontsize='xx-small', legend_fontweight='normal'
			)
	else:
		palette=sns.husl_palette(ncolor)
		sc.pl.umap(x, color=splitby, palette=palette, frameon=False, show=False, save=f'{bname}_umap_{splitby}_wolabel.png')
		sc.pl.umap(x, color=splitby, palette=palette, frameon=False, show=False, save=f'{bname}_umap_{splitby}_ondata.png'
			, legend_loc='on data', legend_fontsize='xx-small', legend_fontweight='normal'
			)

print(dir(x))
print(vars(x))
sc.write(filename=f'{bname}.h5ad', adata=x)
x.obs['barcode']=x.obs.index
x.obs.to_csv(f'{bname}_obs.txt.gz', sep='\t', index=False)
x.var['symbol']=x.var.index
x.var.to_csv(f'{bname}_var.txt.gz', sep='\t', index=False)
