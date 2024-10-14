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

# Default parameters
default_params={
	'scrnah5adsubsetbyvaluecounts': {
		'condaenv': None,
		'label': 'sampleid',
		'ncell': 10,
	},
	'scrnascvih5ad': {
		'condaenv': None,
		'batchkey': 'sampleid',
		'nlayer': 2,
		'nlatent': 30,
		'ntop': 10000,
		'flavor': 'seurat',
		'seed': 12345,
		'epoch': None,
		'gpu': None,
		'normcounts': False,
	},
	'scrnascanpycombinerawcountsscvi': {
		'condaenv': None,
		'obs': ['_scvi_batch', '_scvi_labels', '_scvi_local_l_mean', '_scvi_local_l_var', 'leiden'],
		'invert': True,
	},
	'scrnah5ad2normscale': {
		'condaenv': None,
		'scale': False,
	},
	'scrnah5adumapby': {
		'condaenv': None,
		'width': 5,
		'height': 5,
		'notitle': True,
		'label': [],
	},
}

# Assign default params to config if not existing
for key, paramdict in default_params.items():
	if key not in config:
		config[key]=paramdict
	else:
		for subkey, value in paramdict.items():
			if subkey not in config[key]:
				config[key][subkey]=value
			else:
				print(f"Info: {key} and {subkey} are in config. So, skipping default assignment.")

# Overwrite condaenv
if config['condaenv']=='None':
	config['condaenv']=None
if config['condaenv']:
	for key, paramdict in config.items():
		if isinstance(paramdict, dict) and 'condaenv' in paramdict:
			paramdict['condaenv']=config['condaenv']

# save parameters
with open(f"config_{nowtimestr}.yaml", 'w') as f:
	yaml.dump(config, f, sort_keys=False)

# Build command string for each step
# +
scrnah5adsubsetbyvaluecounts_cmd=[
	f"scrnah5adsubsetbyvaluecounts.sh",
	f"-l {config['scrnah5adsubsetbyvaluecounts']['label']}",
	f"-n {config['scrnah5adsubsetbyvaluecounts']['ncell']}",
	]
if config['scrnah5adsubsetbyvaluecounts']['condaenv']:
	scrnah5adsubsetbyvaluecounts_cmd+=[f"-e {config['scrnah5adsubsetbyvaluecounts']['condaenv']}"]
scrnah5adsubsetbyvaluecounts_cmd=' '.join(scrnah5adsubsetbyvaluecounts_cmd)
# -

# +
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
# -

# +
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
# -

# +
scrnah5ad2normscale_cmd=[
	f"scrnah5ad2normscale.sh",
	]
if config['scrnah5ad2normscale']['condaenv']:
	scrnah5ad2normscale_cmd+=[f"-e {config['scrnah5ad2normscale']['condaenv']}"]
if 'scale' in config['scrnah5ad2normscale'] and config['scrnah5ad2normscale']['scale']:
	scrnah5ad2normscale_cmd+=['-s']
scrnah5ad2normscale_cmd=' '.join(scrnah5ad2normscale_cmd)
# -

# +
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
# -

