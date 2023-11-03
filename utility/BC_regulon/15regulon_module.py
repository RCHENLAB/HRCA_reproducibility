import seaborn as sns
import matplotlib.pyplot as plt

import dill
#work_dir = 'pbmc_tutorial'

celltype="BC"
work_dir=f'/storage/chentemp/wangj/scenic/{celltype}'

scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))

# Correlation between region based regulons and gene based regulons
import pandas as pd
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

#cormt=scplus_obj.uns['eRegulon_AUC']['Gene_based'][keep_gene].corr()
cormt=scplus_obj.uns['eRegulon_AUC']['Region_based'][keep_region].corr()
#cormt.to_csv(f'{work_dir}/Regulon_correlation.txt.gz',sep="\t")

cormt.to_csv(f'{work_dir}/Regulon_correlation_region.txt.gz',sep="\t")
from scipy.spatial import distance
from scipy.cluster import hierarchy
import numpy as np
#correlations = cormt

from scipy.cluster import hierarchy
# Compute the linkage matrix
linkage_matrix = hierarchy.linkage(cormt, method='average')


from scipy.cluster.hierarchy import fcluster
# Cut the dendrogram to create clusters
num_clusters = 10  # Adjust this based on your desired number of modules
clusters = fcluster(linkage_matrix, num_clusters, criterion='maxclust')

# Reorder rows and columns based on cluster assignments
sorted_idx = np.argsort(clusters)
reordered_matrix = cormt.iloc[sorted_idx, sorted_idx]

reordered_matrix.to_csv(f'{work_dir}/Regulon_correlation_reorder_region.txt.gz',sep="\t")

