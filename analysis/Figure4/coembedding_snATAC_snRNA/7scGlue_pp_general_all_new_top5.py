import scvelo as scv
import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import anndata
import sys
#import pytorch-gpu
import torch
import pandas as pd
import numpy as np
import episcanpy.api as epi

#scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (7, 7)

cell=sys.argv[1]
label=sys.argv[2] #"BC_Macaque"
out=sys.argv[3]
label_atac=sys.argv[4]  #"predictedGroup"
label_atac1=sys.argv[5]
in_dir=f'/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1/h5ad/'

sample_list="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1/snATAC_list"
sample=[]
with open(sample_list, "r") as sl:
	for line in sl:
		line1=line.strip()
		sample.append(line1)


adata_list=[]

for sam in sample:
	path=f'{in_dir}/{sam}_flt.h5ad'
	adata=scv.read(path, cache=False)
	adata.obs["sample"]=sam
	adata.obs.index=adata.obs.cellname
	adata_list.append(adata)

atac=anndata.concat(adata_list)

atac.write(f'{out}/atac61_sample.h5ad')

BC_table="/storage/chenlab/Users/junwang/human_meta/data/ATAC_cell_archr_021623"

BC=pd.read_csv(BC_table,sep="\t")

cellname=atac.obs.index.intersection(BC.index)

atac.obs.loc[cellname, label_atac]=BC.loc[cellname, label_atac]


NN_table="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_crossmap_epi/NN_atac_obs.txt.gz"

NN=pd.read_csv(NN_table,sep="\t",index_col=0)

NN_atac=NN[NN.domain=="atac"]

cellname=atac.obs.index.intersection(NN_atac.index)

atac.obs.loc[cellname, label_atac]=NN_atac.loc[cellname, label_atac1]

rna=scv.read("/storage/chenlab/Users/junwang/human_meta/data/ref/major_clean_major_scvi_Cluster_clean.h5ad",cache=False)

rna=rna[rna.obs[label] != "unassigned"]

rna.layers["counts"] = rna.X.copy()

atac.layers["counts"] = atac.X.copy()

rna.X
rna.X.data

sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3", subset=True)

sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")

sc.pp.neighbors(rna, metric="cosine")
sc.tl.umap(rna)
sc.pl.embedding(rna, basis="X_umap", color=[label],ncols=1,frameon=False,save=f'_{cell}_rna.png', palette="tab20")

atac.X
atac.X.data
print(f'start atac hvg')

print(atac.shape)
print(np.max(atac.X))
epi.pp.binarize(atac)
print(np.max(atac.X))

epi.pp.filter_features(atac, min_cells=10)

epi.pp.cal_var(atac)

var_annot = atac.var.sort_values(ascending=False, by ='variability_score')

peak_num=int(atac.shape[1] * 0.1)

high_variable=var_annot[0:peak_num].index

atac.var["highly_variable"]=False

atac.var.loc[high_variable,"highly_variable"]=True

print(atac.shape)
print(f'end atac hvg')

scglue.data.lsi(atac, n_components=100, n_iter=15, use_highly_variable=True) #, random_state=0)
#scglue.data.lsi(atac, n_components=100, n_iter=15)

print(f'start atac neigh')

sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
sc.tl.umap(atac)

print(f'end atac umap')

sc.pl.embedding(atac, basis="X_umap",color=[label_atac], ncols=1,frameon=False,save=f'_{cell}_atac.png', palette="tab20")

rna.var.head()

scglue.data.get_gene_annotation(
    rna, gtf="/storage/singlecell/jinli/resource/cellranger/human/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
    gtf_by="gene_name"
)

rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()

atac.var_names[:5]

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
atac.var.head()

guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
guidance

scglue.graph.check_graph(guidance, [rna, atac])

atac.var.head()

rna.write(f'{out}/{cell}_rna-pp.h5ad', compression="gzip")
atac.write(f'{out}/{cell}_atac-pp.h5ad', compression="gzip")
nx.write_graphml(guidance, f'{out}/{cell}_guidance.graphml.gz')

