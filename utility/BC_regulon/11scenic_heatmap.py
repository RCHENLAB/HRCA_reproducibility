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

work_dir=f'/storage/chentemp/wangj/scenic/{celltype}' # f'/storage/chentemp/wangj/scenic/major'
tmp_dir=f'/storage/chentemp/jw2'

sc.settings.set_figure_params(dpi=200, frameon=False, figsize=(5, 5), facecolor='white')



import dill
#work_dir = 'pbmc_tutorial'
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))

import numpy as np
n_targets = [int(x.split('(')[1].replace('r)', '')) for x in scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
rho = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
adj_pval = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()

thresholds = {
        'rho': [-0.75, 0.70],
        'n_targets': 0
}

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
        figsize = (10, 5),
#        orientation = 'vertical',
        orientation = 'horizontal',
	save=f'{work_dir}/heatmap_dotplot.png')


