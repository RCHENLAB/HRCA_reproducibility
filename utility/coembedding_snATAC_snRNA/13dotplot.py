import scanpy as sc
import matplotlib.pyplot as plt


marker_genes_dict =[
	 'PDE6A',
	 'ARR3',
	 'VSX2',	 
	 'GAD1',
	 'RBPMS',
         'ONECUT2',
         'RLBP1',
       	 'GFAP',
  	 'CD74' ]
 
#	 'VSX2',
#	 'ARR3',
#	 'ONECUT2',
#	 'RLBP1',
#	 'CD74',
#	 'NEFM',
#	 'BEST1',
#	 'PDE6A']

#atac_file="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_downsample2000/major_GeneScoreMatrix.h5ad"
atac_file="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_downsample2000/final_major_downsample2000_GeneScoreMatrix.h5ad"
atac=sc.read_h5ad(atac_file)

#atac
sc.pp.normalize_total(atac, target_sum=1e4)

sc.pp.log1p(atac)
sc.pp.highly_variable_genes(atac, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(atac, max_value=10)
sc.pp.pca(atac) #calculate PCA embedding

#sc.pl.dotplot(atac, marker_genes_dict, groupby='celltype',save=f'major_dot_plots_atac.png', use_raw=False,standard_scale="group")


rna_file="/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_major_downsample2000.h5ad"
rna0=sc.read_h5ad(rna_file)

rna=rna0[rna0.obs.majorclass!="RPE"].copy()

#rna
sc.pp.normalize_total(rna, target_sum=1e4)

sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(rna, max_value=10)
sc.pp.pca(rna) #calculate PCA embedding

#deg=sc.tl.rank_genes_groups(rna)

#sc.pl.dotplot(rna, marker_genes_dict, groupby='majorclass',save=f'major_dot_plots_rna.png',use_raw=False,standard_scale="group")

fig, (ax1, ax2) =plt.subplots(1,2, figsize=(10,4)) # gridspec_kw={'wspace':0.9})
#ax1_dict=sc.pl.dotplot(atac, marker_genes_dict, groupby='celltype',show=False, use_raw=False,standard_scale="group",ax=ax1, categories_order=['FMB', 'OFFx', 'DB1', 'DB2', 'DB3a', 'DB3b', 'IMB', 'GB', 'BB', 'DB4a', 'DB4b', 'DB5', 'DB6', 'RBC'],  cmap="Blues", title="atac", colorbar_title=f'Mean gene score\nin group')

ax1_dict=sc.pl.dotplot(atac, marker_genes_dict, groupby='celltype',show=False, use_raw=False,standard_scale="group",ax=ax1, categories_order=['Rod', 'Cone','BC', 'AC', 'RGC', 'HC', 'MG', 'Astrocyte','Microglia'], cmap="Blues", colorbar_title=f'Mean gene score\nin group',swap_axes=True)
ax2_dict=sc.pl.dotplot(rna, marker_genes_dict, groupby='majorclass',show=False,use_raw=False,standard_scale="group",ax=ax2,  categories_order=['Rod', 'Cone','BC', 'AC', 'RGC', 'HC', 'MG', 'Astrocyte','Microglia'], swap_axes=True)

fig.tight_layout()

#fig.savefig('/storage/chenlab/Users/junwang/human_meta/data/major_downsample2000_dotplot.png', dpi=500)
fig.savefig('/storage/chenlab/Users/junwang/human_meta/data/major_downsample2000_dotplot.pdf', dpi=500)

