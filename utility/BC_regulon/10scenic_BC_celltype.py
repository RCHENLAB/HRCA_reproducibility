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

celltype=sys.argv[1]
model_num=int(sys.argv[2])
label=sys.argv[3]
metacell_num=int(sys.argv[4])
#work_dir = f'/storage/chenlab/Users/junwang/human_meta/data/scenic/{celltype}'
#tmp_dir = f'/storage/chenlab/Users/junwang/'

#celltype="major"
#model_num=32


work_dir=f'/storage/chentemp/wangj/scenic/{celltype}' # f'/storage/chentemp/wangj/scenic/major'
tmp_dir=f'/storage/chentemp/jw2'

sc.settings.set_figure_params(dpi=200, frameon=False, figsize=(5, 5), facecolor='white')

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



#rna_h5ad="/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_major_downsample2000.h5ad"   #f'/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_{celltype}_downsample1000.h5ad'
rna_h5ad=f'/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_{celltype}_downsample1000.h5ad'
rna=sc.read_h5ad(rna_h5ad)
#label="cluster2"
rna.obs["celltype"]=rna.obs[label]

if celltype == "BC":
	rna.obs["celltype"]=rna.obs["celltype"].replace(BC_ct)

if celltype == "Cone":
	rna.obs["celltype"]=rna.obs["celltype"].replace(Cone_ct)

ncell=rna.n_obs*0.005

sc.pp.filter_genes(rna, min_cells=ncell)


#######
#create cistopic

#atac_h5ad=f'/storage/chenlab/Users/junwang/human_meta/data/ref/snATAC_{celltype}_downsample1000.h5ad'
atac_h5ad=f'/storage/chenlab/Users/junwang/human_meta/data/ref/snATAC_{celltype}_downsample1000.h5ad'
#atac_h5ad=f'/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_downsample2000/final_major_downsample2000_PeakMatrix.h5ad'

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

#cistopic_obj = create_cistopic_object(fragment_matrix=sparse_X, split_pattern="", cell_names=atac.obs_names.values, region_names=atac.var_names.str.replace("-",":",1).values, project="" )
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
#exit()




######inferring candidate enhancer
from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)
#######
###calculate DARs per cell type
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='celltype', var_features=variable_regions)

#####save the result
if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))
import pickle
pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'wb'))

import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')

db_fpath="/storage/chenlab/Users/junwang/reference"

rankings_db = os.path.join(db_fpath, 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_fpath, 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(db_fpath, 'motifs-v10-nr.hgnc-m0.00001-o0.0.tbl')


#we will run pycistarget using the run_pycistarget wrapper function

if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))

from scenicplus.wrappers.run_pycistarget import run_pycistarget
run_pycistarget(
    region_sets = region_sets,
    species = 'homo_sapiens',
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 7,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = 'v10nr_clust',
    )

import dill
menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))

menr['DEM_topics_otsu_All'].DEM_results('Topic2')

##########infer enhancer-driven Gene Regulatory Networks
adata=rna
from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata, #adata.raw.to_adata(),
    cisTopic_obj = cistopic_obj,
    ACC_prefix = 'ACC_',
    GEX_prefix = 'GEX_',
    multi_ome_mode=False,
    nr_metacells=int(metacell_num),
    menr = menr,
    key_to_group_by='celltype',
    bc_transform_func = lambda x: f'{x}'  #-10x_pbmc' #function to convert scATAC-seq barcodes to scRNA-seq ones
)
#scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
scplus_obj

if not os.path.exists(os.path.join(work_dir, 'scenicplus')):
    os.makedirs(os.path.join(work_dir, 'scenicplus'))

dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)

#####select ensembl id
ensembl_version_dict = {'105': 'http://www.ensembl.org',
                        '104': 'http://may2021.archive.ensembl.org/',
                        '103': 'http://feb2021.archive.ensembl.org/',
                        '102': 'http://nov2020.archive.ensembl.org/',
                        '101': 'http://aug2020.archive.ensembl.org/',
                        '100': 'http://apr2020.archive.ensembl.org/',
                        '99': 'http://jan2020.archive.ensembl.org/',
                        '98': 'http://sep2019.archive.ensembl.org/',
                        '97': 'http://jul2019.archive.ensembl.org/',
                        '96': 'http://apr2019.archive.ensembl.org/',
                        '95': 'http://jan2019.archive.ensembl.org/',
                        '94': 'http://oct2018.archive.ensembl.org/',
                        '93': 'http://jul2018.archive.ensembl.org/',
                        '92': 'http://apr2018.archive.ensembl.org/',
                        '91': 'http://dec2017.archive.ensembl.org/',
                        '90': 'http://aug2017.archive.ensembl.org/',
                        '89': 'http://may2017.archive.ensembl.org/',
                        '88': 'http://mar2017.archive.ensembl.org/',
                        '87': 'http://dec2016.archive.ensembl.org/',
                        '86': 'http://oct2016.archive.ensembl.org/',
                        '80': 'http://may2015.archive.ensembl.org/',
                        '77': 'http://oct2014.archive.ensembl.org/',
                        '75': 'http://feb2014.archive.ensembl.org/',
                        '54': 'http://may2009.archive.ensembl.org/'}

import pybiomart as pbm
def test_ensembl_host(scplus_obj, host, species):
    dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
    print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
    return ov

n_overlap = {}
for version in ensembl_version_dict.keys():
    print(f'host: {version}')
    try:
        n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'hsapiens')
    except:
        print('Host not reachable')
v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")

biomart_host = "http://sep2019.archive.ensembl.org/"


import scenicplus
from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *

import dill
#work_dir = 'pbmc_tutorial'
#scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))

filter_genes(scplus_obj, min_pct = 0.5)
filter_regions(scplus_obj, min_pct = 0.5)

from scenicplus.wrappers.run_scenicplus import run_scenicplus
try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['celltype'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = '/storage/chenlab/Users/junwang/human_meta/data/scenic/utoronto_human_tfs_v_1.01.txt',
        save_path = os.path.join(work_dir, 'scenicplus'),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = False,
        export_to_UCSC_file = False,
        path_bedToBigBed = '/storage/chen/home/jw29/software',
        n_cpu = 7,
        _temp_dir = os.path.join(tmp_dir, 'ray_spill'))
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
    raise(e)


########
for attr in dir(scplus_obj.uns['eRegulons'][0]):
    if not attr.startswith('_'):
        print(f"{attr}: {getattr(scplus_obj.uns['eRegulons'][0], attr) if not type(getattr(scplus_obj.uns['eRegulons'][0], attr)) == list else getattr(scplus_obj.uns['eRegulons'][0], attr)[0:5]}")

scplus_obj.uns['eRegulon_metadata'].head()

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

import dill
#work_dir = 'pbmc_tutorial'
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))

from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
apply_std_filtering_to_eRegulons(scplus_obj)

scplus_obj.uns['eRegulon_metadata_filtered'].head()

from scenicplus.eregulon_enrichment import score_eRegulons
region_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus/region_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
gene_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus/gene_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
score_eRegulons(scplus_obj,
                ranking = region_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type= 'region',
                auc_threshold = 0.05,
                normalize = False,
                n_cpu = 7)
score_eRegulons(scplus_obj,
                gene_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type = 'gene',
                auc_threshold = 0.05,
                normalize= False,
                n_cpu = 7)

from scenicplus.dimensionality_reduction import run_eRegulons_tsne, run_eRegulons_umap
run_eRegulons_umap(
    scplus_obj = scplus_obj,
    auc_key = 'eRegulon_AUC_filtered',
    reduction_name = 'eRegulons_UMAP', #overwrite previously calculated UMAP
)


run_eRegulons_tsne(
    scplus_obj = scplus_obj,
    auc_key = 'eRegulon_AUC_filtered',
    reduction_name = 'eRegulons_tSNE', #overwrite previously calculated tSNE
)

#######
from scenicplus.dimensionality_reduction import plot_metadata_given_ax
import matplotlib.pyplot as plt
import seaborn as sns
#%matplotlib inline

#specify color_dictionary

#color_dict = {
#    'Rod': "#E9842C",
#    'Cone': "#F8766D",
#    'BC': "#BC9D00",
#    'AC': "#00C0B4",
#    'Astrocyte': "#9CA700",
#    'MG': "#6FB000" ,
#    'Microglia': "#00B813",
#    'HC': "#00BD61",
#    'RGC': '#00C08E'
#}


color_dict = {
    'DB3b': "#E9842C",
    'DB4b': "#F8766D",
    'DB4a': "#BC9D00",
    'DB5': "#00C0B4",
    'DB2': "#9CA700",
    'IMB': "#6FB000" ,
    'FMB': "#00B813",
    'RBC': "#00BD61",
    'OFFx': '#00C08E',
    'DB6': '#00BDD4',
    'DB1': '#00A7FF',
    'BB': '#7F96FF',
    'GB': '#E26EF7',
    'DB3a': '#FF62BF'}



fig, axs = plt.subplots(ncols=2, figsize = (16, 8))
plot_metadata_given_ax(
    scplus_obj=scplus_obj,
    ax = axs[0],
    reduction_name = 'eRegulons_UMAP',
    variable = 'celltype', #note the GEX_ prefix, this metadata originated from the gene expression metadata (on which we did the cell type annotation before)
    color_dictionary={'celltype': color_dict}
)
plot_metadata_given_ax(
    scplus_obj=scplus_obj,
    ax = axs[1],
    reduction_name = 'eRegulons_tSNE',
    variable = 'celltype', #note the GEX_ prefix, this metadata originated from the gene expression metadata (on which we did the cell type annotation before)
    color_dictionary={'celltype': color_dict}
)
fig.tight_layout()
sns.despine(ax = axs[0]) #remove top and right edge of axis border
sns.despine(ax = axs[1]) #remove top and right edge of axis border
#plt.show()
plt.savefig(f'{work_dir}/eRegulons.png')

#########
from scenicplus.cistromes import TF_cistrome_correlation, generate_pseudobulks

generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Gene_based')


generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Region_based')

TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'celltype',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Gene_based',
            out_key = 'filtered_gene_based')
TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'celltype',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Region_based',
            out_key = 'filtered_region_based')

scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].head()

import numpy as np
n_targets = [int(x.split('(')[1].replace('r)', '')) for x in scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
rho = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
adj_pval = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()

thresholds = {
        'rho': [-0.75, 0.70],
        'n_targets': 0
}
import seaborn as sns
fig, ax = plt.subplots(figsize = (10, 5))
sc = ax.scatter(rho, n_targets, c = -np.log10(adj_pval), s = 5)
ax.set_xlabel('Correlation coefficient')
ax.set_ylabel('nr. target regions')
#ax.hlines(y = thresholds['n_targets'], xmin = min(rho), xmax = max(rho), color = 'black', ls = 'dashed', lw = 1)
ax.vlines(x = thresholds['rho'], ymin = 0, ymax = max(n_targets), color = 'black', ls = 'dashed', lw = 1)
ax.text(x = thresholds['rho'][0], y = max(n_targets), s = str(thresholds['rho'][0]))
ax.text(x = thresholds['rho'][1], y = max(n_targets), s = str(thresholds['rho'][1]))
sns.despine(ax = ax)
fig.colorbar(sc, label = '-log10(adjusted_pvalue)', ax = ax)
plt.savefig(f'{work_dir}/cistrome_correlation.png')

selected_cistromes = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].loc[
        np.logical_or(
                scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] > thresholds['rho'][1],
                scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] < thresholds['rho'][0]
        )]['Cistrome'].to_list()


selected_eRegulons = [x.split('_(')[0] for x in selected_cistromes]
selected_eRegulons_gene_sig = [
        x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()
        if x.split('_(')[0] in selected_eRegulons]
selected_eRegulons_region_sig = [
        x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()
        if x.split('_(')[0] in selected_eRegulons]
#save the results in the scenicplus object
scplus_obj.uns['selected_eRegulon'] = {'Gene_based': selected_eRegulons_gene_sig, 'Region_based': selected_eRegulons_region_sig}
print(f'selected: {len(selected_eRegulons_gene_sig)} eRegulons')


dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)

from scenicplus.plotting.dotplot import heatmap_dotplot
heatmap_dotplot(
        scplus_obj = scplus_obj,
        size_matrix = scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'], #specify what to plot as dot sizes, target region enrichment in this case
        color_matrix = scplus_obj.to_df('EXP'), #specify  what to plot as colors, TF expression in this case
        scale_size_matrix = True,
        scale_color_matrix = True,
        group_variable = 'celltype',
        subset_eRegulons = scplus_obj.uns['selected_eRegulon']['Gene_based'],
        index_order = ['DB3b', 'DB4b', 'DB4a', 'DB5', 'DB2', 'IMB', 'FMB', 'RBC', 'OFFx', 'DB6', 'DB1','BB','GB','DB3a'],
        figsize = (5, 20),
        orientation = 'vertical',
	save=f'{work_dir}/heatmap_dotplot.png')

from scenicplus.RSS import *
regulon_specificity_scores(
        scplus_obj,
        variable = 'celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_keys = ['Region_based'],
        selected_regulons = [x for x in scplus_obj.uns['selected_eRegulon']['Region_based'] if '-' not in x],
        out_key_suffix = '_filtered')

plot_rss(scplus_obj, 'celltype_filtered', num_columns=4, top_n=10, figsize = (10, 10), save=f'{work_dir}/rss.png')

######add######
from scenicplus.RSS import *
regulon_specificity_scores(
        scplus_obj,
        variable = 'celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_keys = ['Gene_based'],
        selected_regulons = [x for x in scplus_obj.uns['selected_eRegulon']['Gene_based'] if '-' not in x],
        out_key_suffix = '_filtered')

plot_rss(scplus_obj, 'celltype_filtered', num_columns=4, top_n=10, figsize = (10, 10), save=f'{work_dir}/rss_genebased.png')


#############3

flat_list = lambda t: [item for sublist in t for item in sublist]
selected_markers = list(set(flat_list(
    [scplus_obj.uns['RSS']['celltype_filtered'].loc[celltype].sort_values(ascending = False).head(10).index.to_list()
    for celltype in scplus_obj.uns['RSS']['celltype_filtered'].index])))

from scenicplus.plotting.correlation_plot import *


region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Region_based',
        use_plotly = False,
        selected_regulons = selected_markers,
        signature_key = 'eRegulon_signatures_filtered',
        figsize = (10, 10), return_data = True, vmax = 0.5, cmap = 'plasma',save=f'{work_dir}/jaccard_heatmap.png')

from scenicplus.plotting.correlation_plot import *
correlation_heatmap(scplus_obj,
                    auc_key = 'eRegulon_AUC',
                    signature_keys = ['Gene_based'],
                    selected_regulons = scplus_obj.uns['selected_eRegulon']['Gene_based'],
                    fcluster_threshold = 0.1,
                    fontsize = 3, save=f'{work_dir}/correlation_heatmap.png')

