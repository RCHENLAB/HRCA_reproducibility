import sys
import scvelo as scv
import anndata
import pandas as pd
import os
import numpy as np

datafile="/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_other_all2022"
#batch   donor   file
#Chen_a_10x3_Lobe_19_D003_NeuN   Chen_19_D003    /storage/singlecell/jinli/wkfl/metaanalysis/prvdata/lobe/scrnaseurath5ad2h5seurat/10x3_Lobe_19_D003_NeuN.h5ad
donor=sys.argv[1] #"/storage/chenlab/Users/junwang/human_meta/data/donor30"
dir1="/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_raw_batch_new_clean"

#dir1="/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_raw_batch_new"
#os.mkdir(dir1)
ref_dir="/storage/chenlab/Users/junwang/human_meta/data/ref"
cell=["AC","BC","HC","NN","PR","RGC"]
cell_list=[]
for c in cell:
	fn=f'{ref_dir}/{c}.obs.gz'
	cl=pd.read_csv(fn,sep="\t")
	cell_list=cell_list+cl.iloc[:,0].tolist()

def rename_gene(adata):
	for i in range(0, adata.var.features.shape[0]):
		if adata.var.features[i] in gene_dict:
			adata.var.features[i] = gene_dict[adata.var.features[i]]
	return adata

def read_h5ad_list(file_list,donor):
	adata_list=[]
	with open (file_list,"r") as fl:
		for line in fl:
			info=line.strip().split()
			if(info[0] == donor):
				adata=scv.read(info[-1])
				adata.var_names=adata.var.index
				adata.layers["raw"] = adata.X
				adata1=adata[adata.obs.index.isin(cell_list)]
				adata_list.append(adata1)
	adata_full=anndata.concat(adata_list,join="inner") 
	return adata_full

def grouped_obs_cpm(adata, group_key, layer=None, gene_symbols=None):
	if layer is not None:
		getX= lambda x: x.layers[layer]
	else:
		getX = lambda x: x.raw.X
	if gene_symbols is not None:
		new_idx = adata.var[idx]
	else:
		new_idx = adata.var_names
	grouped = adata.obs.groupby(group_key)
	out = pd.DataFrame(
		np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
		columns=list(grouped.groups.keys()),
		index=adata.var_names
	)
	for group, idx in grouped.indices.items():
		X = getX(adata[idx])
		total_count=X.sum(dtype=np.float64)
#		factor=1000000.0
#		out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))/total_count*factor
		out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))
	return(out)

with open(donor, "r") as da:
	for line in da:
		da1=line.strip()
		adata_full = read_h5ad_list(datafile,da1)
		genexp=grouped_obs_cpm(adata=adata_full, group_key="scpred_prediction", layer="raw")
		fileout=f'{dir1}/{da1}.txt.gz'
		genexp.to_csv(fileout,sep="\t")
	
