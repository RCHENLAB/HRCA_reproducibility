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


work_dir=f'/storage/chentemp/wangj/scenic/{celltype}' # f'/storage/chentemp/wangj/scenic/major'
tmp_dir=f'/storage/chentemp/jw_{celltype}'

sc.settings.set_figure_params(dpi=500, frameon=False, figsize=(5, 5), facecolor='white')

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

import dill
#work_dir = 'pbmc_tutorial'
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))

# Correlation between region based regulons and gene based regulons
import pandas
df1 = scplus_obj.uns['eRegulon_AUC']['Gene_based'].copy()
df2 = scplus_obj.uns['eRegulon_AUC']['Region_based'].copy()
df1.columns = [x.split('_(')[0] for x in df1.columns]
df2.columns = [x.split('_(')[0] for x in df2.columns]
correlations = df1.corrwith(df2, axis = 0)
correlations = correlations[abs(correlations) > 0.6]
# Kepp only R2G +
keep = [x for x in correlations.index if '+_+' in x] + [x for x in correlations.index if '-_+' in x]
# Keep extended if not direct
extended = [x for x in keep if 'extended' in x]
direct = [x for x in keep if not 'extended' in x]
keep_extended = [x for x in extended if not x.replace('extended_', '') in direct]
keep = direct + keep_extended
# Keep regulons with more than 10 genes
keep_gene = [x for x in scplus_obj.uns['eRegulon_AUC']['Gene_based'].columns if x.split('_(')[0] in keep]
keep_gene = [x for x in keep_gene if (int(x.split('_(')[1].replace('g)', '')) > 10)]
keep_all = [x.split('_(')[0] for x in keep_gene]
keep_region = [x for x in scplus_obj.uns['eRegulon_AUC']['Region_based'].columns if x.split('_(')[0] in keep]
scplus_obj.uns['selected_eRegulons'] = {}
scplus_obj.uns['selected_eRegulons']['Gene_based'] = keep_gene
scplus_obj.uns['selected_eRegulons']['Region_based'] = keep_region

len(keep_gene)

from scenicplus.plotting.correlation_plot import *
correlation_heatmap(scplus_obj,
                    auc_key = 'eRegulon_AUC',
                    signature_keys = ['Gene_based'],
                    selected_regulons = scplus_obj.uns['selected_eRegulons']['Gene_based'],
                    fcluster_threshold = 0.1,
                    fontsize = 7, save=f'{work_dir}/correlation_heatmap.pdf', figsize = (20, 20))


from scenicplus.plotting.correlation_plot import *
jaccard_heatmap(scplus_obj,
                    gene_or_region_based = 'Gene_based',
                    signature_key = 'eRegulon_signatures',
                    selected_regulons = scplus_obj.uns['selected_eRegulons']['Gene_based'],
                    fcluster_threshold = 0.1,
                    fontsize = 7,
                    method='intersect', save=f'{work_dir}/correlation_heatmap_eRegulon_gene.pdf', figsize = (20, 20), cmap = 'plasma')

from scenicplus.plotting.correlation_plot import *
jaccard_heatmap(scplus_obj,
                    gene_or_region_based = 'Region_based',
                    signature_key = 'eRegulon_signatures',
                    selected_regulons = scplus_obj.uns['selected_eRegulons']['Region_based'],
                    fcluster_threshold = 0.1,
                    fontsize = 7,
                    method='intersect', save=f'{work_dir}/correlation_heatmap_eRegulon_region.pdf', figsize = (20, 20),  cmap = 'plasma')


region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Region_based',
#        use_plotly = False,
        fontsize = 10,
        fcluster_threshold = 0.1,
        selected_regulons = scplus_obj.uns['selected_eRegulon']['Region_based'],
        signature_key = 'eRegulon_signatures_filtered',
#        plot_dendrogram =True,
        figsize = (12, 10), return_data = True, vmax = 0.5, cmap = 'plasma',save=f'{work_dir}/jaccard_heatmap_selected_Regulon_region.pdf')



region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Gene_based',
#        use_plotly = False,
        fontsize = 10,
        fcluster_threshold = 0.1,
        selected_regulons = scplus_obj.uns['selected_eRegulon']['Gene_based'],
        signature_key = 'eRegulon_signatures_filtered',
        figsize = (12, 10), return_data = True, vmax = 0.5, cmap = 'magma',save=f'{work_dir}/jaccard_heatmap_selected_Regulon_gene.pdf')


