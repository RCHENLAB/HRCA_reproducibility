import yaml
import subprocess
from pathlib import Path

# Parameters
infile, outdir, bname, nowtimestr=config['infile'].split(','), config['outdir'], config['bname'].split(','), config['nowtimestr']
config['infile'], config['bname']=infile, bname
files=dict(zip(bname, infile))

# for rules using a function
def get_file(wildcards):
	return files[wildcards.name]

if config['condaenv']=='None':
	config['condaenv']=None

# Default parameters
if 'scrnah5adduplicateobs' not in config:
	config['scrnah5adduplicateobs']={
		'condaenv': config['condaenv'],
		'saturn_label': 'celltype',
		'batch_label': 'sampleid',
		'global': False,
	}
if config['condaenv']: # set a higher priority for the option
	config['scrnah5adduplicateobs']['condaenv']=config['condaenv']

scrnah5adduplicateobs_cmd=[
	f"scrnah5adduplicateobs",
	]
if config['scrnah5adduplicateobs']['condaenv']:
	scrnah5adduplicateobs_cmd+=[f"-e {config['scrnah5adduplicateobs']['condaenv']}"]
if config['scrnah5adduplicateobs']['global']: # using global label for all species
	scrnah5adduplicateobs_cmd+=[
		f"-s {config['scrnah5adduplicateobs']['saturn_label']}",
		f"-t saturn_label",
		f"-s {config['scrnah5adduplicateobs']['batch_label']}",
		f"-t batch_label",
	]
scrnah5adduplicateobs_cmd=' '.join(scrnah5adduplicateobs_cmd)



if 'scrnah5adsubsetsamplingbykey' not in config:
	config['scrnah5adsubsetsamplingbykey']={
		'skip': False,
		'condaenv': config['condaenv'],
		'key': 'saturn_label',
		'nsample': 2000,
		'seed': 12345,
	}
if config['condaenv']: # set a higher priority for the option
	config['scrnah5adsubsetsamplingbykey']['condaenv']=config['condaenv']

scrnah5adsubsetsamplingbykey_cmd=[
	f"scrnah5adsubsetsamplingbykey.sh",
	f"-k {config['scrnah5adsubsetsamplingbykey']['key']}",
	f"-n {config['scrnah5adsubsetsamplingbykey']['nsample']}",
	f"-s {config['scrnah5adsubsetsamplingbykey']['seed']}",
	]
if config['scrnah5adsubsetsamplingbykey']['condaenv']:
	scrnah5adsubsetsamplingbykey_cmd+=[f"-e {config['scrnah5adsubsetsamplingbykey']['condaenv']}"]
scrnah5adsubsetsamplingbykey_cmd=' '.join(scrnah5adsubsetsamplingbykey_cmd)



if 'scrnah5adsubsetbyvaluecounts' not in config:
	config['scrnah5adsubsetbyvaluecounts']={
		'skip': False,
		'condaenv': config['condaenv'],
		'label': 'batch_label',
		'ncell': 10,
	}
if config['condaenv']: # set a higher priority for the option
	config['scrnah5adsubsetbyvaluecounts']['condaenv']=config['condaenv']

scrnah5adsubsetbyvaluecounts_cmd=[
	f"scrnah5adsubsetbyvaluecounts.sh",
	f"-l {config['scrnah5adsubsetbyvaluecounts']['label']}",
	f"-n {config['scrnah5adsubsetbyvaluecounts']['ncell']}",
	]
if config['scrnah5adsubsetbyvaluecounts']['condaenv']:
	scrnah5adsubsetbyvaluecounts_cmd+=[f"-e {config['scrnah5adsubsetbyvaluecounts']['condaenv']}"]
scrnah5adsubsetbyvaluecounts_cmd=' '.join(scrnah5adsubsetbyvaluecounts_cmd)




if 'saturntraincrossspecies' not in config:
	config['saturntraincrossspecies']={
		'label': 'saturn_label',
		'batchkey': 'batch_label',
		'condaenv': config['condaenv'],
		'ngene': 5000,
		'nhvg': 8000,
		'hvgspan': 1.0,
		'mapfile': None,
		'epoch': 50,
		'gpu': 0,
		'seed': 12345,
	}
if config['condaenv']: # set a higher priority for the option
	config['saturntraincrossspecies']['condaenv']=config['condaenv']

saturntraincrossspecies_cmd=[
	f"saturntraincrossspecies",
	f"-l {config['saturntraincrossspecies']['label']}",
	f"-k {config['saturntraincrossspecies']['batchkey']}",
	f"-m {config['saturntraincrossspecies']['ngene']}",
	f"-n {config['saturntraincrossspecies']['nhvg']}",
	f"-v {config['saturntraincrossspecies']['hvgspan']}",
	f"-p {config['saturntraincrossspecies']['epoch']}",
	f"-g {config['saturntraincrossspecies']['gpu']}",
	f"-s {config['saturntraincrossspecies']['seed']}",
	]
if config['saturntraincrossspecies']['condaenv']:
	saturntraincrossspecies_cmd+=[f"-e {config['saturntraincrossspecies']['condaenv']}"]
if config['saturntraincrossspecies']['mapfile']:
	saturntraincrossspecies_cmd+=[f"-f {config['saturntraincrossspecies']['mapfile']}"]
saturntraincrossspecies_cmd=' '.join(saturntraincrossspecies_cmd)


if 'saturntrainh5ad2umap' not in config:
	config['saturntrainh5ad2umap']={
		'condaenv': config['condaenv'],
		'label': ['batch_labels', 'labels', 'labels2', 'ref_labels', 'species'],
		'seed': 12345,
		'resolution': 0.5,
		'neighbor': 15,
		'npc': None,
	}
if config['condaenv']: # set a higher priority for the option
	config['saturntrainh5ad2umap']['condaenv']=config['condaenv']

saturntrainh5ad2umap_cmd=[
	f"saturntrainh5ad2umap",
	f"-s {config['saturntrainh5ad2umap']['seed']}",
	f"-r {config['saturntrainh5ad2umap']['resolution']}",
	f"-n {config['saturntrainh5ad2umap']['neighbor']}",
	]
if config['saturntrainh5ad2umap']['condaenv']:
	saturntrainh5ad2umap_cmd+=[f"-e {config['saturntrainh5ad2umap']['condaenv']}"]
if config['saturntrainh5ad2umap']['npc']:
	saturntrainh5ad2umap_cmd+=[f"-p"]
if 'label' in config['saturntrainh5ad2umap'] and len(config['saturntrainh5ad2umap']['label'])>0:
	saturntrainh5ad2umap_cmd+=[f"-l {s}" for s in config['saturntrainh5ad2umap']['label']]
saturntrainh5ad2umap_cmd=' '.join(saturntrainh5ad2umap_cmd)



if 'saturnh5adumap2seuratdimplotsplitby' not in config:
	config['saturnh5adumap2seuratdimplotsplitby']={
		'condaenv': config['condaenv'],
		'group': ['labels', 'labels2', 'ref_labels'],
		'split': 'species',
		'height': 6,
		'width': 6,
		'legendwidth': 4,
	}
if config['condaenv']:
	config['saturnh5adumap2seuratdimplotsplitby']['condaenv']=config['condaenv']

saturnh5adumap2seuratdimplotsplitby_cmd=[
	f"saturnh5adumap2seuratdimplotsplitby",
	f"-s {config['saturnh5adumap2seuratdimplotsplitby']['split']}",
	f"-H {config['saturnh5adumap2seuratdimplotsplitby']['height']}",
	f"-W {config['saturnh5adumap2seuratdimplotsplitby']['width']}",
	f"-L {config['saturnh5adumap2seuratdimplotsplitby']['legendwidth']}",
	]
if config['saturnh5adumap2seuratdimplotsplitby']['condaenv']:
	saturnh5adumap2seuratdimplotsplitby_cmd+=[f"-e {config['saturnh5adumap2seuratdimplotsplitby']['condaenv']}"]
if 'group' in config['saturnh5adumap2seuratdimplotsplitby'] and len(config['saturnh5adumap2seuratdimplotsplitby']['group'])>0:
	saturnh5adumap2seuratdimplotsplitby_cmd+=[f"-g {g}" for g in config['saturnh5adumap2seuratdimplotsplitby']['group']]
saturnh5adumap2seuratdimplotsplitby_cmd=' '.join(saturnh5adumap2seuratdimplotsplitby_cmd)



if 'saturnh5adumap2seuratclustree' not in config:
	config['saturnh5adumap2seuratclustree']={
		'condaenv': config['condaenv'],
		'group': 'labels',
		'color': 'species',
		'reorder': 'F',
		'numericorder': 'F',
		'height': 7.5,
		'width': 6,
	}
if config['condaenv']:
	config['saturnh5adumap2seuratclustree']['condaenv']=config['condaenv']

saturnh5adumap2seuratclustree_cmd=[
	f"saturnh5adumap2seuratclustree",
	f"-g {config['saturnh5adumap2seuratclustree']['group']}",
	f"-c {config['saturnh5adumap2seuratclustree']['color']}",
	f"-r {config['saturnh5adumap2seuratclustree']['reorder']}",
	f"-n {config['saturnh5adumap2seuratclustree']['numericorder']}",
	f"-H {config['saturnh5adumap2seuratclustree']['height']}",
	f"-W {config['saturnh5adumap2seuratclustree']['width']}",
	]
if config['saturnh5adumap2seuratclustree']['condaenv']:
	saturnh5adumap2seuratclustree_cmd+=[f"-e {config['saturnh5adumap2seuratclustree']['condaenv']}"]
saturnh5adumap2seuratclustree_cmd=' '.join(saturnh5adumap2seuratclustree_cmd)


if 'scrnah5adsplitby' not in config:
	config['scrnah5adsplitby']={
		'condaenv': config['condaenv'],
		'splitby': 'species',
	}
if config['condaenv']: # set a higher priority for the option
	config['scrnah5adsplitby']['condaenv']=config['condaenv']

scrnah5adsplitby_cmd=[
	f"scrnah5adsplitby.sh",
	f"-s {config['scrnah5adsplitby']['splitby']}",
	]
if config['scrnah5adsplitby']['condaenv']:
	scrnah5adsplitby_cmd+=[f"-e {config['scrnah5adsplitby']['condaenv']}"]
scrnah5adsplitby_cmd=' '.join(scrnah5adsplitby_cmd)



if 'saturnh5adlatent2cntclassifier' not in config:
	config['saturnh5adlatent2cntclassifier']={
		'condaenv': config['condaenv'],
		# 'model': ['LogisticRegression', 'SVC', 'KNeighborsClassifier', 'RandomForestClassifier', 'SVClinear', 'SGDClassifier', 'GradientBoostingClassifier', 'MultinomialNB', 'NeighborhoodComponentsAnalysis'],
		# 'model': ['LogisticRegression', 'SVC', 'KNeighborsClassifier', 'RandomForestClassifier', 'SVClinear', 'SGDClassifier', 'MultinomialNB'],
		'model': ['LogisticRegression', 'KNeighborsClassifier', 'RandomForestClassifier'],
		'trainlabel': 'labels',
		'predictlabel': 'labels',
		'threshold': [0.55, 0.9],
		'seed': 12345,
	}
if config['condaenv']: # set a higher priority for the option
	config['saturnh5adlatent2cntclassifier']['condaenv']=config['condaenv']

saturnh5adlatent2cntclassifier_train_cmd=[
	f"saturnh5adlatent2cntclassifier train",
	f"-c {config['saturnh5adlatent2cntclassifier']['trainlabel']}",
	f"-s {config['saturnh5adlatent2cntclassifier']['seed']}",
	]
if config['saturnh5adlatent2cntclassifier']['condaenv']:
	saturnh5adlatent2cntclassifier_train_cmd+=[f"-e {config['saturnh5adlatent2cntclassifier']['condaenv']}"]
saturnh5adlatent2cntclassifier_train_cmd=' '.join(saturnh5adlatent2cntclassifier_train_cmd)

saturnh5adlatent2cntclassifier_predict_cmd=[
	f"saturnh5adlatent2cntclassifier predict",
	f"-l {config['saturnh5adlatent2cntclassifier']['predictlabel']}",
	]
if config['saturnh5adlatent2cntclassifier']['condaenv']:
	saturnh5adlatent2cntclassifier_predict_cmd+=[f"-e {config['saturnh5adlatent2cntclassifier']['condaenv']}"]
saturnh5adlatent2cntclassifier_predict_cmd=' '.join(saturnh5adlatent2cntclassifier_predict_cmd)

saturnh5adlatent2cntclassifier_predictproba_cmd=[
	f"saturnh5adlatent2cntclassifier predictproba",
	f"-l {config['saturnh5adlatent2cntclassifier']['predictlabel']}",
	]
if config['saturnh5adlatent2cntclassifier']['condaenv']:
	saturnh5adlatent2cntclassifier_predictproba_cmd+=[f"-e {config['saturnh5adlatent2cntclassifier']['condaenv']}"]
saturnh5adlatent2cntclassifier_predictproba_cmd=' '.join(saturnh5adlatent2cntclassifier_predictproba_cmd)


if 'rtable2sankeydiagram' not in config:
	config['rtable2sankeydiagram']={
		'condaenv': config['condaenv'],
		'iteration': 32,
	}
if config['condaenv']: # set a higher priority for the option
	config['rtable2sankeydiagram']['condaenv']=config['condaenv']

rtable2sankeydiagram_cmd=[
	f"rtable2sankeydiagram",
	f"-t {config['rtable2sankeydiagram']['iteration']}",
	]
if config['rtable2sankeydiagram']['condaenv']:
	rtable2sankeydiagram_cmd+=[f"-e {config['rtable2sankeydiagram']['condaenv']}"]
rtable2sankeydiagram_cmd=' '.join(rtable2sankeydiagram_cmd)


if 'tsvfileaddfromcolumn' not in config:
	config['tsvfileaddfromcolumn']={
		'key': ['predict', 'predict_max'],
	}
tsvfileaddfromcolumn_cmd=[
	f"tsvfileaddfromcolumn",
	]
if 'key' in config['tsvfileaddfromcolumn'] and len(config['tsvfileaddfromcolumn']['key'])>0:
	tsvfileaddfromcolumn_cmd+=[f"-k {k}" for k in config['tsvfileaddfromcolumn']['key']]
tsvfileaddfromcolumn_cmd=' '.join(tsvfileaddfromcolumn_cmd)



if 'saturncntclassifierpred2tab' not in config:
	config['saturncntclassifierpred2tab']={
		'xlabel': ['saturn_label', 'group1', 'group2', 'group3', 'cluster2'],
		'ylabel': ['predict', 'predict_max'],
	}
saturncntclassifierpred2tab_cmd=[
	f"saturncntclassifierpred2tab",
	]
if 'xlabel' in config['saturncntclassifierpred2tab'] and len(config['saturncntclassifierpred2tab']['xlabel'])>0:
	saturncntclassifierpred2tab_cmd+=[f"-x {k}" for k in config['saturncntclassifierpred2tab']['xlabel']]
if 'ylabel' in config['saturncntclassifierpred2tab'] and len(config['saturncntclassifierpred2tab']['ylabel'])>0:
	saturncntclassifierpred2tab_cmd+=[f"-y {k}" for k in config['saturncntclassifierpred2tab']['ylabel']]
saturncntclassifierpred2tab_cmd=' '.join(saturncntclassifierpred2tab_cmd)



# debug parameters
with open(f"config_{nowtimestr}.yaml", 'w') as f:
	yaml.dump(config, f, sort_keys=False)
