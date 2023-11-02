from itertools import chain
import matplotlib.pyplot as plt
import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams
import sys
import numpy as np
import random

random.seed(10)

indir=sys.argv[1]
celltype=sys.argv[2]
rna_batch=sys.argv[3] #sampleid
atac_batch=sys.argv[4] #sample
label=sys.argv[5]
atac_label=sys.argv[6]
outdir=sys.argv[7]
#scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (7, 7)

rna = ad.read_h5ad(f'{indir}/{celltype}_rna-pp.h5ad')
rna.obs["domain"]="rna"
rna.obs["cell_type"]=rna.obs[label]
rna.obs[atac_batch]=rna.obs[rna_batch]
atac = ad.read_h5ad(f'{indir}/{celltype}_atac-pp.h5ad')
atac.obs["domain"]="atac"
atac.obs["cell_type"]=atac.obs[atac_label]
guidance = nx.read_graphml(f'{indir}/{celltype}_guidance.graphml.gz')

scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca",
    use_batch=rna_batch,
    use_cell_type="cell_type"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi",
    use_batch=atac_batch,
    use_cell_type="cell_type"
)

guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()


glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory": "glue"} #,init_kws={"random_seed": 0}   
)

glue.save(f'{outdir}/{celltype}_glue.dill')

#glue.save(f'{indir}/{celltype}_glue.dill')
# glue = scglue.models.load_model("glue.dill")

dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, guidance_hvf
)
dx

sns_plot=sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")
#plt.savefig(f'{indir}/{celltype}_integration_consistency.png')
plt.savefig(f'{outdir}/{celltype}_integration_consistency.png')

rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

combined = ad.concat([rna, atac])

sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.leiden(combined)
sc.tl.umap(combined)
#sc.pl.umap(combined, color=["scpred_prediction", "domain"], wspace=0.65)
#sc.pl.embedding(combined, basis="X_umap", color=[label,atac_label,"cell_type", "domain"],ncols=1,frameon=False,save=f'_{celltype}_combined.png', palette="tab20")
#sc.pl.embedding(combined, basis="X_umap", color=["leiden","cell_type", "domain"],ncols=1,frameon=False,save=f'_{celltype}_combined.png', palette="tab20")
#sc.pl.embedding(combined, basis="X_umap", color=["leiden","cell_type", "domain",atac_batch],ncols=2,frameon=False,save=f'_{celltype}_combined_cpu.png', palette="tab20")
#sc.pl.embedding(combined, basis="X_umap", color=["cell_type","leiden", "domain",atac_batch],ncols=2,frameon=False,save=f'_{celltype}_combined.png', palette="tab20")
sc.pl.embedding(combined, basis="X_umap", color=["cell_type","leiden", atac_batch],ncols=2,frameon=False,save=f'_{celltype}_combined.png', palette="tab20")

sc.pl.embedding(combined, basis="X_umap", color=["leiden"],ncols=1,legend_loc="on data",frameon=False,save=f'_{celltype}_combined_leiden.png', palette="tab20")

sc.pl.embedding(combined, basis="X_umap", color=["domain"],ncols=1,frameon=False,save=f'_{celltype}_combined_domain.png',palette=["lightblue","blue"])

sc.pl.embedding(combined, basis="X_umap", color=["cell_type"],frameon=False,save=f'_{celltype}_combined.png', palette="tab20")

sc.pl.embedding(combined, basis="X_umap", color=[atac_batch],frameon=False,save=f'_{celltype}_combined_sample.png', palette="tab20")


df=combined[(combined.obs.domain=="rna")].obs.groupby(["leiden","cell_type"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat.to_csv(f'{outdir}/{celltype}_conf_mat_rna.csv',header=True, index=True)


df=combined[(combined.obs.domain=="atac")].obs.groupby(["leiden","cell_type"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat.to_csv(f'{outdir}/{celltype}_conf_mat_atac.csv',header=True, index=True)


df = combined.obs.groupby(["leiden", "cell_type"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]

fig=plt.figure(figsize=(10, 15))
fig=plt.pcolor(conf_mat)
fig= plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
fig= plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
fig=plt.xlabel("Predicted")
fig=plt.ylabel("Observed")
#plt.savefig(f'{indir}/{celltype}_conf_mat.png',format="png",facecolor="w")
#conf_mat.to_csv(f'{indir}/{celltype}_conf_mat_table.csv',header=True, index=True)

plt.savefig(f'{outdir}/{celltype}_conf_mat.png',format="png",facecolor="w")
conf_mat.to_csv(f'{outdir}/{celltype}_conf_mat_table.csv',header=True, index=True)

feature_embeddings = glue.encode_graph(guidance_hvf)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :5]

rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()

#combined.write(f'{indir}/{celltype}_combined-emb.h5ad', compression="gzip")

#rna.write(f'{indir}/{celltype}_rna-emb.h5ad', compression="gzip")
#atac.write(f'{indir}/{celltype}_atac-emb.h5ad', compression="gzip")
#nx.write_graphml(guidance_hvf, f'{indir}/{celltype}_guidance-hvf.graphml.gz')

combined.write(f'{outdir}/{celltype}_combined-emb.h5ad', compression="gzip")

rna.write(f'{outdir}/{celltype}_rna-emb.h5ad', compression="gzip")
atac.write(f'{outdir}/{celltype}_atac-emb.h5ad', compression="gzip")
nx.write_graphml(guidance_hvf, f'{outdir}/{celltype}_guidance-hvf.graphml.gz')


