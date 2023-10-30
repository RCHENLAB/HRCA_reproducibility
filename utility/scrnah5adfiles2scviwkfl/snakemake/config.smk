from snakemake.utils import min_version
min_version("7.0")
import yaml

# Parameters
infile, outdir, bname, nowtimestr=config['infile'].split(','), config['outdir'], config['bname'].split(','), config['nowtimestr']
config['infile'], config['bname']=infile, bname
files=dict(zip(bname, infile))

# global wildcards constraints for all rules
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards
wildcard_constraints:
  name="[^/]+"

# for rules using a function
def get_file(wildcards):
	return files[wildcards.name]

if config['condaenv']=='None':
	config['condaenv']=None

# Default parameters

if 'scrnah5adsubsetbyvaluecounts' not in config:
	config['scrnah5adsubsetbyvaluecounts']={
		'condaenv': config['condaenv'],
		'label': 'sampleid',
		'ncell': 10,
	}
if config['condaenv']:
	config['scrnah5adsubsetbyvaluecounts']['condaenv']=config['condaenv']

scrnah5adsubsetbyvaluecounts_cmd=[
	f"scrnah5adsubsetbyvaluecounts.sh",
	f"-l {config['scrnah5adsubsetbyvaluecounts']['label']}",
	f"-n {config['scrnah5adsubsetbyvaluecounts']['ncell']}",
	]
if config['scrnah5adsubsetbyvaluecounts']['condaenv']:
	scrnah5adsubsetbyvaluecounts_cmd+=[f"-e {config['scrnah5adsubsetbyvaluecounts']['condaenv']}"]
scrnah5adsubsetbyvaluecounts_cmd=' '.join(scrnah5adsubsetbyvaluecounts_cmd)


if 'scrnascvih5ad' not in config:
	config['scrnascvih5ad']={
		'condaenv': config['condaenv'],
		'batchkey': 'sampleid',
		'nlayer': 2,
		'nlatent': 30,
		'ntop': 10000,
		'flavor': 'seurat',
		'seed': 12345,
		'epoch': None,
		'gpu': None,
		'normcounts': False,
	}
if config['condaenv']:
	config['scrnascvih5ad']['condaenv']=config['condaenv']

scrnascvih5ad_cmd=[
	f"scrnascvih5ad.sh",
	f"-k {config['scrnascvih5ad']['batchkey']}",
	f"-l {config['scrnascvih5ad']['nlayer']}",
	f"-t {config['scrnascvih5ad']['nlatent']}",
	f"-n {config['scrnascvih5ad']['ntop']}",
	f"-f {config['scrnascvih5ad']['flavor']}",
	f"-s {config['scrnascvih5ad']['seed']}",
	f"-p {config['scrnascvih5ad']['epoch']}",
	]
if config['scrnascvih5ad']['condaenv']:
	scrnascvih5ad_cmd+=[f"-e {config['scrnascvih5ad']['condaenv']}"]
if 'gpu' in config['scrnascvih5ad'] and config['scrnascvih5ad']['gpu']:
	scrnascvih5ad_cmd+=[f"-g {config['scrnascvih5ad']['gpu']}"]
if 'normcounts' in config['scrnascvih5ad'] and config['scrnascvih5ad']['normcounts']:
	scrnascvih5ad_cmd+=['--normcounts']
scrnascvih5ad_cmd=' '.join(scrnascvih5ad_cmd)

# scrnascanpycombinerawcountsscvi
if 'scrnascanpycombinerawcountsscvi' not in config:
	config['scrnascanpycombinerawcountsscvi']={
		'condaenv': config['condaenv'],
		'obs': ['_scvi_batch', '_scvi_labels', '_scvi_local_l_mean', '_scvi_local_l_var'],
		'invert': True,
	}
if config['condaenv']:
	config['scrnascanpycombinerawcountsscvi']['condaenv']=config['condaenv']

scrnascanpycombinerawcountsscvi_cmd=[
	f"scrnascanpycombinerawcountsscvi.sh",
	]
if config['scrnascanpycombinerawcountsscvi']['condaenv']:
	scrnascanpycombinerawcountsscvi_cmd+=[f"-e {config['scrnascanpycombinerawcountsscvi']['condaenv']}"]
if 'obs' in config['scrnascanpycombinerawcountsscvi'] and len(config['scrnascanpycombinerawcountsscvi']['obs'])>0:
	scrnascanpycombinerawcountsscvi_cmd+=['-s '+s for s in config['scrnascanpycombinerawcountsscvi']['obs']]
if 'invert' in config['scrnascanpycombinerawcountsscvi'] and config['scrnascanpycombinerawcountsscvi']['invert']:
	scrnascanpycombinerawcountsscvi_cmd+=['-n']
scrnascanpycombinerawcountsscvi_cmd=' '.join(scrnascanpycombinerawcountsscvi_cmd)



# scrnah5ad2normscale
if 'scrnah5ad2normscale' not in config:
	config['scrnah5ad2normscale']={
		'condaenv': config['condaenv'],
		'scale': False,
	}
if config['condaenv']:
	config['scrnah5ad2normscale']['condaenv']=config['condaenv']

scrnah5ad2normscale_cmd=[
	f"scrnah5ad2normscale.sh",
	]
if config['scrnah5ad2normscale']['condaenv']:
	scrnah5ad2normscale_cmd+=[f"-e {config['scrnah5ad2normscale']['condaenv']}"]
if 'scale' in config['scrnah5ad2normscale'] and config['scrnah5ad2normscale']['scale']:
	scrnah5ad2normscale_cmd+=['-s']
scrnah5ad2normscale_cmd=' '.join(scrnah5ad2normscale_cmd)



# scrnah5adumapby
if 'scrnah5adumapby' not in config:
	config['scrnah5adumapby']={
		'condaenv': config['condaenv'],
		'width': 5,
		'height': 5,
		'notitle': True,
		'label': ['sampleid', 'leiden', 'nCount_RNA', 'nFeature_RNA', 'percent.mt'],
	}
if config['condaenv']:
	config['scrnah5adumapby']['condaenv']=config['condaenv']

scrnah5adumapby_cmd=[
	f"scrnah5adumapby.sh",
	f"-W {config['scrnah5adumapby']['width']}",
	f"-H {config['scrnah5adumapby']['height']}",
	]
if config['scrnah5adumapby']['condaenv']:
	scrnah5adumapby_cmd+=[f"-e {config['scrnah5adumapby']['condaenv']}"]
if config['scrnah5adumapby']['notitle']:
	scrnah5adumapby_cmd+=[f"-t"]
if 'label' in config['scrnah5adumapby'] and len(config['scrnah5adumapby']['label'])>0:
	scrnah5adumapby_cmd+=['-l '+l for l in config['scrnah5adumapby']['label']]
scrnah5adumapby_cmd=' '.join(scrnah5adumapby_cmd)

# debug parameters
with open(f"config_{nowtimestr}.yaml", 'w') as f:
	yaml.dump(config, f, sort_keys=False)
