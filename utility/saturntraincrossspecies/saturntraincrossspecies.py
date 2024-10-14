#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
import pandas as pd
from exeBash.exeBash import exeBash
from pathlib import Path
import click
import random
import socket

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-s', '--seed', type=click.INT, default=12345, show_default=True, help='Random seed.')
@click.option('-l', '--label', type=click.STRING, default='celltype', help='Label of input data.')
@click.option('-g', '--gpu', type=click.INT, default=0, show_default=True, help='A GPU device.')
@click.option('-m', '--ngene', type=click.INT, default=2000, show_default=True, help='The number of macrogenes.')
@click.option('-n', '--nhvg', type=click.INT, default=8000, show_default=True, help='The number of HVGs.')
@click.option('-v', '--hvgspan', type=click.FLOAT, default=1.0, show_default=True, help='HVG span for seurat_v3.')
@click.option('-f', '--mapfile', type=click.Path(exists=True, resolve_path=True), help='Cell type map file.')
@click.option('-p', '--epoch', type=click.INT, default=50, show_default=True, help='Number of epoch.')
@click.option('-k', '--batchkey', type=click.STRING, help='Non-speceis batch key. E.g., sampleid.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, seed, label, gpu, ngene, nhvg, hvgspan, mapfile, epoch, batchkey, infile):
	"""
Run SATURN train for cross-species analysis.

\b
INFILE is a sample TSV file with three columns: species<TAB>path<TAB>embedding_path. The header names are fixed.
1. `species`: species
2. `path`: .h5ad file path
3. `embedding_path`: protein embeddings

\b
```E.g.,
$ cat fnameinfo.txt
species	path	embedding_path
human	/storage/singlecell/jinli/mwe/SATURNmwe/example/BC/download/scrnah5adsubsetsamplingbykey/scrnah5adduplicateobs/BC.h5ad	/storage/singlecell/jinli/resource/saturn/protein_embeddings/protein_embeddings_export/ESM2/human_embedding.torch
mouse	/storage/singlecell/jinli/mwe/SATURNmwe/example/BC/download/scrnah5adsubsetsamplingbykey/scrnah5adduplicateobs/mm_concat47_inner.h5ad	/storage/singlecell/jinli/resource/saturn/protein_embeddings/protein_embeddings_export/ESM2/mouse_embedding.torch
```

\b
Example:
  slurmtaco.sh -n g01 -m 120G -- saturntraincrossspecies -d "$outdir" -b BC_human_mouse -e saturn -l saturn_label -k sampleid -- fnameinfo.txt

\b
Note:
  1. `-f|-mapfile` is used to add scoring metrics when training when available.
  2. Be aware that barcode is lost for species, however, the order of cells are retained. So, it can recover other metadata by the order.

\b
See also:
  Depends:
    Github/SATURN

\b
Date: 2023/02/07
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	Path(outdir).mkdir(parents=True, exist_ok=True)
	config=pd.read_csv(infile, sep='\t', header=0)
	config['path']=[Path().absolute() / p for p in config['path']]
	config.to_csv(f"{outdir}/{bname}.csv", index=False)

	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	os.chdir(outdir)
	script=f"{absdir}/SATURN/train-saturn.py"
	exprs=[
		f"{script}",
		f"--in_data '{bname}.csv'",
		f"--device 'cuda'",
		f"--work_dir '.'",
		f"--seed {seed}",
		f"--in_label_col '{label}'",
		f"--ref_label_col '{label}'",
		f"--num_macrogenes {ngene}",
		f"--hv_genes {nhvg}",
		f"--hv_span {hvgspan}",
		f"--epochs {epoch}",
		f"--centroids_init_path '{bname}_centroids.pkl'",
		]

	if gpu>0:
		exprs+=[f"--device_num {gpu}"]

	if mapfile:
		exprs+=[
			f"--score_adata",
			f"--ct_map_path {mapfile}",
			]

	if batchkey:
		exprs+=[
			f"--non_species_batch_col '{batchkey}'",
			]

	return exeBash.callback(exprs, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
