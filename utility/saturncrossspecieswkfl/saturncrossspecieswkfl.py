#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
import shutil
from exeBash.exeBash import exeBash
from pathlib import Path
import datetime
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment, e.g., sccaf.')
@click.option('-c', '--configfile', type=click.Path(exists=False, resolve_path=True), help='Configuration file in the YAML format.')
@click.option('-t', '--numthreads', type=click.INT, default=4, show_default=True, help='Number of threads.')
@click.option('-b', '--bname', multiple=True, type=click.STRING, help='Bnames for infile. Default: basenames of infiles.')
@click.option('-n', '--dryrun', is_flag=True, help='Dry-run.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True), required=True, nargs=-1)
def main(outdir, condaenv, configfile, numthreads, bname, dryrun, infile):
	"""
To perform cross-species analysis by SATURN.

INFILE are configuration YAML file[s]. E.g.,

\b
```BC_mouse_macaque.yaml
- species: human
  path: /storage/singlecell/jinli/wkfl/atlashumanprj/integration/snRNA/cross_map/preproc/scrnah5adsplitby/snRNA_BC.h5ad
  embedding_path: /storage/singlecell/jinli/resource/saturn/protein_embeddings/protein_embeddings_export/ESM2/human_embedding.torch
  saturn_label: cluster2
  batch_label: sampleid
  classifier: predict
- species: mouse
  path: /storage/singlecell/jinli/wkfl/atlashumanprj/integration/snRNA/cross_species/preproc/mouse/clean/BC.h5ad
  embedding_path: /storage/singlecell/jinli/resource/saturn/protein_embeddings/protein_embeddings_export/ESM2/mouse_embedding.torch
  saturn_label: subclass
  batch_label: sampleid
  classifier: train
- species: macaque
  path: /storage/singlecell/jinli/wkfl/atlashumanprj/integration/snRNA/cross_species/preproc/macaque/scrnah5adsplitby/full_BC.h5ad
  embedding_path: /storage/singlecell/jinli/resource/saturn/protein_embeddings/protein_embeddings_export/ESM2/macaca_fascicularis_embedding.torch
  saturn_label: group1
  batch_label: replicateid
  classifier: train
```

\b
Example:
  slurmtaco.sh -n g01 -t 4 -m 60G -- saturncrossspecieswkfl -d "$outdir" -e saturn -t 4 -n -- "$indir"/*.yaml
  # slurmtaco.sh -n g01 -t 4 -m 60G -- saturncrossspecieswkfl -d "$outdir" -e saturn -t 4 -c config.yaml -- "$indir"/*.yaml

\b
Note:
  1. `-e|--condaenv` will set the conda environment for all internal steps, and it has a higher priority over the `config.yaml`.

\b
See also:
  Steps:
    scrnah5adduplicateobs
    scrnah5adsubsetsamplingbykey.sh
    scrnah5adsubsetbyvaluecounts.sh
    saturntraincrossspecies
    saturntrainh5ad2umap
    saturnh5adumap2seuratdimplotsplitby
    saturnh5adumap2seuratclustree
    scrnah5adsplitby
    saturnh5adlatent2cntclassifier
    rtable2sankeydiagram
    html2png: run on the local machine
    tsvfileaddfromcolumn
    tsvfileaddfromcolumn/saturncntclassifierpred2tab
    tsvfileaddfromcolumn/rtable2sankeydiagram

\b
Date: 2023/02/10
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if len(bname)!=len(infile):
		bname=[Path(file).stem for file in infile]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	nowtimestr=datetime.datetime.now().strftime('%y%m%d_%H%M%S')

	for file in [
			f"{absdir}/snakemake/Snakefile",
			f"{absdir}/snakemake/config.smk",
			]:
		shutil.copy2(file, outdir) # copy Snakfile and config.smk

	# common command string
	cmdstr=[
		f"snakemake",
		f"-C outdir='{outdir}' infile='{','.join(infile)}' bname='{','.join(bname)}' condaenv='{condaenv}' nowtimestr='{nowtimestr}'",
		f"-d {outdir}",
		f"-j {numthreads}",
		]
	if configfile:
		cmdstr+=[f"--configfile {configfile}"]

	if dryrun: # dry-run
		cmdstr+=[f"--dry-run --printshellcmds"]

	else: # running
		cmdstr+=[
			f"--printshellcmds --debug-dag --skip-script-cleanup --verbose",
			]

	os.chdir(outdir)
	return exeBash.callback(cmdstr=cmdstr, condaenv=None, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
