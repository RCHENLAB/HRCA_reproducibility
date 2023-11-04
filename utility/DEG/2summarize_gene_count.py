import sys
import scvelo as scv
import anndata
import pandas as pd
import os
import numpy as np

datafile="/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_other_all2023_mac_lobe_batch"
#batch   donor   file
#Chen_a_10x3_Lobe_19_D003_NeuN   Chen_19_D003    /storage/singlecell/jinli/wkfl/metaanalysis/prvdata/lobe/scrnaseurath5ad2h5seurat/10x3_Lobe_19_D003_NeuN.h5ad
donor="/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_other_all2023_mac_lobe_batch_donor_rmChang"
#os.mkdir(dir1)

cell_list=[]
ref_dir="/storage/singlecell/jinli/wkfl/atlashumanprj/integration/snRNA"
cell=["cellbrowser/v1/scrnah5adaddmetadatamerge/snRNA","cross_map/preproc/scrnah5adrenameclusterfromobs/scrnah5adfiles2scviwkfl/scRNA"]
for c in cell:
	fn=f'{ref_dir}/{c}_obs.txt.gz'
	cl=pd.read_csv(fn,sep="\t")
	cell_list=cell_list+cl.iloc[:,-1].tolist()


def read_h5ad_list(file_list,donor):
	adata_list=[]
	print(donor)
	with open (file_list,"r") as fl:
		for line in fl:
			info=line.strip().split()
			if(info[2] == donor):
				adata=scv.read(info[-1])
				if (info[1].find("chen_d") == -1 ) and (info[1].find("chen_e") == -1 ) and (info[1].find("chen_f") == -1 ) :
					adata1=adata[adata.obs.index.isin(cell_list)]
				else :
					adata1=adata
				adata_list.append(adata1)
	adata_full=anndata.concat(adata_list,join="inner") 
	return adata_full

celltype=["AC","BC","Cone","HC","MG","RGC","Rod","Astrocyte","Microglia","RPE"]
cellnum={}

list1="/storage/chenlab/Users/junwang/human_meta/data/donor_all_reform_batch_new_snRNA_clean2023_all"

out=open(list1,"w")
with open(donor, "r") as da:
	for line in da:
		da1=line.strip()
		data = read_h5ad_list(datafile,da1)
		for cell in celltype:
			if cell in data.obs.scpred_prediction.value_counts():
				if da1 not in cellnum:
					cellnum[da1]={}
				cellnum[da1][cell] = data.obs.scpred_prediction.value_counts()[cell]
				num_tmp=data.obs.scpred_prediction.value_counts()[cell]
				out.write(f'{da1}\t{cell}\t{num_tmp}\n')		
			else:
				if da1 not in cellnum:
					cellnum[da1]={}
				cellnum[da1][cell] = 0
				num_tmp=0
				out.write(f'{da1}\t{cell}\t{num_tmp}\n')

out.close()
df=pd.DataFrame.from_dict(cellnum)

list="/storage/chenlab/Users/junwang/human_meta/data/donor_all_batch_new_snRNA_clean2023_all"

df.to_csv(f'{list}_celltype_num')
	
