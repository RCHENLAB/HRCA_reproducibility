#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exeR.exeR import exeR
from pathlib import Path
import random
import socket
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-k', '--key', type=click.STRING, required=True, help='Key for groups.')
@click.option('-1', '--group1', type=click.STRING, required=True, help='Group 1.')
@click.option('-2', '--group2', type=click.STRING, required=True, help='Group 2. This is the baseline.')
@click.option('-l', '--label', type=click.STRING, required=True, help='Key for a cluster variable. E.g., majorclass.')
@click.option('-r', '--rowmean', type=click.FLOAT, default=1.0, show_default=True, help='Row means to filter out low-count genes.')
@click.option('-c', '--colmean', type=click.FLOAT, default=1.0, show_default=True, help='Column means to filter out low-count samples.')
@click.option('-W', '--width', type=click.FLOAT, default=8, show_default=True, help='Width.')
@click.option('-H', '--height', type=click.FLOAT, default=5, show_default=True, help='Height.')
@click.option('-s', '--breaks', type=click.INT, default=100, show_default=True, help='Breaks.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, key, group1, group2, label, rowmean, colmean, width, height, breaks, infile):
	"""
Perform differential gene expression analysis between two conditions for cell types from a psuedobulk .h5ad file using DESeq2. x.X is raw counts.

INFILE is a .h5ad pseudobulk file.

\b
Example:
  indir=$(mrrdir.sh ../scrnah5adrawcounts2bulkby)
  outdir=$(mrrdir.sh)
  function cmd {
  local f=$1
  local bname=$(basename "$f" .h5ad)
  if fileexists.sh "$f"
  then
    slurmtaco.sh -t 2 -m 20G -- scrnah5ad2bulkdegcondtypebydeseq2 -e deseq2 -d "$outdir" -b "$bname" -k dataset -1 snRNA -2 scRNA -l majorclass -r 100 -c 10 -H 5 -W 8 -s 100 -- "$f"
  fi
  }
  source env_parallel.bash
  env_parallel cmd ::: "$indir"/*.h5ad

\b
See also:
  Upstream:
    scrnah5adrawcounts2bulkby
    scrnah5adX2rowmeanshist
  Depends:
    R/rhdf5
    R/DESeq2

\b
Date: 2023/05/09
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=[
		f"{absdir}/R/{scriptname}.R",
		]
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"key='{key}'",
		f"group1='{group1}'",
		f"group2='{group2}'",
		f"label='{label}'",
		f"rowmean={rowmean}",
		f"colmean={colmean}",
		f"width={width}",
		f"height={height}",
		f"breaks={breaks}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	# os.chdir(outdir)
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
