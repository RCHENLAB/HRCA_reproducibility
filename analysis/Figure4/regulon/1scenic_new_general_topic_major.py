import sys
import os
import scanpy as sc
import numpy as np
import scanpy as sc
import pandas as pd
import pickle
import pycisTopic
from pycisTopic.cistopic_class import *
from scipy import sparse
import dill
import pyranges


celltype="major"
model_num=32


work_dir="/storage/chentemp/wangj/scenic/major" # f'/storage/chentemp/wangj/scenic/major'
tmp_dir=f'/storage/chentemp/jw1'

sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(5, 5), facecolor='white')

BC_ct={
"HBC14": "BB",
"HBC9": "GB",
"HBC6":"DB1",
"HBC4":	"DB2",
"HBC12": "DB3a",
"HBC8": "DB3b",
"HBC7": "DB4a",
"HBC13": "DB4b",
"HBC5": "DB5",
"HBC10": "DB6",
"HBC1":	"FMB",
"HBC3":	"IMB",
"HBC11": "OFFx",
"HBC2":	"RBC"}

Cone_ct={
"Cone1": "ML_Cone",
"Cone2": "ML_Cone",
"Cone3": "S_Cone"}



rna_h5ad="/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_major_downsample2000.h5ad"   #f'/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_{celltype}_downsample1000.h5ad'
rna=sc.read_h5ad(rna_h5ad)
label="majorclass"
rna.obs["celltype"]=rna.obs[label]

if celltype == "BC":
	rna.obs["celltype"]=rna.obs["celltype"].replace(BC_ct)

if celltype == "Cone":
	rna.obs["celltype"]=rna.obs["celltype"].replace(Cone_ct)

ncell=rna.n_obs*0.005

sc.pp.filter_genes(rna, min_cells=ncell)


#######
#create cistopic

atac_h5ad=f'/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_downsample2000/final_major_downsample2000_PeakMatrix.h5ad'

atac=sc.read_h5ad(atac_h5ad)

if celltype == "BC":
	atac.obs["celltype"]=atac.obs["celltype"].replace(BC_ct)

if celltype == "Cone":
	atac.obs["celltype"]=atac.obs["celltype"].replace(Cone_ct)


ncell=atac.n_obs*0.005

sc.pp.filter_genes(atac, min_cells=ncell)

blacklist="/storage/chenlab/Users/junwang/reference/hg38-blacklist.v2.bed.gz"

sparse_X = sparse.csr_matrix(atac.X.T)

cistopic_obj = create_cistopic_object(fragment_matrix=sparse_X, split_pattern="", cell_names=atac.obs_names.values, region_names=atac.var_names.values, project="" )

# Adding cell information
atac_cell_data_meta=atac.obs #atac_cell_data.loc[atac.obs.index,:].copy()

####create cisTopic object
cistopic_obj.add_cell_data(atac_cell_data_meta)
print(cistopic_obj)
##### save object
if not os.path.exists(os.path.join(work_dir, 'scATAC')):
	os.makedirs(os.path.join(work_dir, 'scATAC'))

pickle.dump(cistopic_obj,open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))
##########

########generate topic modeling
models=run_cgs_models(cistopic_obj, n_topics=[2,4,10,16,32,48], n_cpu=5,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None,
                    _temp_dir = os.path.join(tmp_dir + 'ray_spill'))

if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
    os.makedirs(os.path.join(work_dir, 'scATAC/models'))

pickle.dump(models,
            open(os.path.join(work_dir, 'scATAC/models/10x_pbmc_models_500_iter_LDA.pkl'), 'wb'))

from pycisTopic.lda_models import *
model = evaluate_models(models,
                       select_model=model_num,
                       return_model=True,
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False, save=f'{work_dir}/model_evalution.pdf')
cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

##########
from pycisTopic.clust_vis import *
run_umap(cistopic_obj, target  = 'cell', scale=True)
plot_metadata(cistopic_obj, reduction_name = 'UMAP', variables = ['celltype', "Sample"], target='cell', num_columns=4,
                 text_size=10,
                 dot_size=5,
                 figsize=(10,5), save=f'{work_dir}/cisTopic_umap.pdf')
plot_topic(cistopic_obj, reduction_name = 'UMAP', num_columns = 4, save=f'{work_dir}/cisTopic_umap_prob.pdf')
########



