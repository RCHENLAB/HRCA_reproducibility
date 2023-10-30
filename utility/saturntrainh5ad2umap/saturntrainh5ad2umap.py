#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exePython.exePython import exePython
from pathlib import Path
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-l', '--label', type=click.STRING, multiple=True, default=['batch_labels', 'labels', 'labels2', 'ref_labels', 'species'], show_default=True, help='Labels of UMAP.')
@click.option('-s', '--seed', type=click.INT, default=12345, show_default=True, help='Random seed.')
@click.option('-r', '--resolution', type=click.FLOAT, default=0.5, show_default=True, help='Resolution for leiden clustering.')
@click.option('-n', '--neighbor', type=click.INT, default=15, show_default=True, help='The number of local neighbors.')
@click.option('-p', '--npc', type=click.INT, help='The number of PCs to calculate neighbors. Full PCs by default.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, label, seed, resolution, neighbor, npc, infile):
	"""
Plot UMAP of cross-species by `saturntraincrossspecies`.

INFILE is a .h5ad file.

\b
Example:
  indir=$(mrrdir.sh ../saturn_results)
  outdir=$(mrrdir.sh)
  function cmd {
  local f=$1
  local bname=$(basename "$f" .h5ad)
  if fileexists.sh "$f"
  then
  	slurmtaco.sh -t 2 -m 20G -- saturntrainh5ad2umap -d "$outdir" -b "$bname" -e saturn -- "$f"
  fi
  }
  source env_parallel.bash
  env_parallel cmd ::: "$indir"/*.h5ad

\b
Note:
  1. Be sure to use anndata=0.8.0, which is required for SATURN.

\b
See also:
  Upstream:
    saturntraincrossspecies
  Related:
    scrnascanpytlumap.sh
    scrnascanpytlumapnohvg.sh

\b
Date: 2023/02/07
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/python/{scriptname}.py'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"label={label}",
		f"seed={seed}",
		f"resolution={resolution}",
		f"neighbor={neighbor}",
		f"npc={npc}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
