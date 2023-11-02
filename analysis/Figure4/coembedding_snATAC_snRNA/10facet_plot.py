import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

adata=sc.read_h5ad("/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_atac1perc_lr/major_clean_atac1perc_lr_combined-emb.h5ad")

dbt_table="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_atac1perc_lr/emb_dbt.txt.gz"

dbt=pd.read_csv(dbt_table,sep="\t",index_col=0)

cone_table="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_crossmap_epi/Cone/Cone_emb_dbt.txt.gz"

dbt_cone=pd.read_csv(cone_table,sep="\t",index_col=0)

dbt_full=dbt.index.append(dbt_cone.index)

celltype="major_clean"
cells=adata.obs.index.difference(dbt_full)
adata_new=adata[cells,:].copy()
label="cell_type"
adata=adata_new[adata_new.obs[label]!="RPE",:]
# Get unique conditions
conditions = adata.obs['domain'].unique()

# Create figure with subplots for each condition
fig, axes = plt.subplots(1, len(conditions), figsize=(9, 3.5), sharex=True, sharey=True)

# Loop over conditions and plot t-SNE with CD11b color for each
for i, condition in enumerate(conditions):
    # Subset data for current condition
	adata_sub = adata[adata.obs['domain'] == condition]
#	sc.pl.umap(adata_sub, color=label, size=10,  title=condition, show=False, ax=axes[i], legend_loc='on data')
	sc.pl.embedding(adata_sub, basis="X_umap", color=[label], size=10, title="", frameon=False, show=False, ax=axes[i] ) #  ,legend_loc='on data') #, palette="tab20")

#	sc.pl.embedding(adata_sub, basis="X_umap", color=[label], size=10,  title=condition, frameon=False, show=False, ax=axes[i] ) #  ,legend_loc='on data') #, palette="tab20")
   
    # Compute t-SNE
#    sc.tl.tsne(adata_sub)
    
    # Plot t-SNE with CD11b color
    
# Save and show the figure
fig.tight_layout()
fig.savefig('/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_atac1perc_lr/umap_plots.pdf', dpi=500)
#plt.show()
