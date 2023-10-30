#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exeR.exeR import exeR
from pathlib import Path
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-g', '--group', type=click.STRING, default='labels', show_default=True, help='Group for clusters.')
@click.option('-H', '--height', type=click.FLOAT, default=5, show_default=True, help='Height.')
@click.option('-W', '--width', type=click.FLOAT, default=5, show_default=True, help='Width.')
@click.option('-n', '--nfeature', type=click.INT, default=2000, show_default=True, help='Number of top variable features.')
@click.option('-a', '--assay', type=click.Choice(['counts', 'data', 'scale.data']), is_flag=False, flag_value='data', default='scale.data', show_default=True, help='Data assay.')
@click.option('-m', '--method', type=click.Choice(['complete', 'average']), is_flag=False, flag_value='complete', default='average', show_default=True, help='Agglomeration method.')
@click.option('-u', '--uhang', type=click.FLOAT, default=0.1, show_default=True, help='Hang for hclust.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, group, height, width, nfeature, assay, method, uhang, infile):
	"""
Run Seurat::BuildClusterTree() from a .h5ad file. It is assumed for .h5ad by anndata=0.8.0, e.g., SATURN.

INFILE is a .h5ad file.

\b
Example:
  indir=/storage/singlecell/jinli/wkfl/atlashumanprj/integration/snRNA/cross_species/saturncrossspecieswkfl/BC/s2k/scrnah5adsubsetbyvaluecounts/BC_mouse_macaque
  outdir=$(mrrdir.sh)
  function cmd {
  local species=$1
  local group=$2
  local f=$indir/$species.h5ad
  local bname=${species}_${group}
  if fileexists.sh "$f"
  then
  	slurmtaco.sh -t 2 -m 20G -- scrnah5adexpr2hclust -e h5ad08 -d "$outdir" -b "$bname" -g "$group" -- "$f"
  fi
  }
  source env_parallel.bash
  env_parallel --colsep='\\t' cmd <<EOF
  human	cluster2
  mouse	subclass
  macaque	group1
  EOF

\b
Note:
  1. `-a scale.data` works perfectly than `counts` or `data`.

\b
See also:
  Related:
    saturnh5adumap2seuratclustree
  Depends:
    R/rhdf5
    R/Seurat

\b
Date: 2023/03/21
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"group='{group}'",
		f"height={height}",
		f"width={width}",
		f"nfeature={nfeature}",
		f"assay='{assay}'",
		f"method='{method}'",
		f"uhang={uhang}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	# os.chdir(outdir)
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
