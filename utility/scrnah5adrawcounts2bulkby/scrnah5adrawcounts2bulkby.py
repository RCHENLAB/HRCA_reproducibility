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
@click.option('-k', '--key', type=click.STRING, multiple=True, help='Keys to summarize pseudo-bulk.')
@click.option('-H', '--idheader', type=click.STRING, default='bulkid', show_default=True, help='Pseudo-bulk ID header.')
@click.option('-s', '--separator', type=click.STRING, default='::', show_default=True, help='Separator to join keys for ID header.')
@click.option('-f', '--filterlabel', type=click.STRING, help='Label to filter cells. E.g., sampleid. No filtering by default.')
@click.option('-c', '--mincell', type=click.INT, default=1000, show_default=True, help='Minimum cells of `-f|--filterlabel`.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, key, idheader, separator, filterlabel, mincell, infile):
	"""
To calculate pseudo-bulk for samples summarized by keys (`-k|--key`).

INFILE is a .h5ad file.

\b
Example:
  f=/storage/singlecell/jinli/wkfl/atlashumanprj/integration/snRNA/cross_map/preproc/scrnah5adsplitby/snRNA_BC.h5ad
  bname=$(basename "$f" .h5ad)
  outdir=$(mrrdir.sh)
  slurmtaco.sh -t 2 -m 20G -- scrnah5adrawcounts2bulkby -d "$outdir" -b "$bname" -k sampleid -k cluster2 -H bulkid -f sampleid -c 1000 -- "$f"

\b
Note:
  1. Raw counds will be aggregated per group.
  2. The resulting count matrix might be good for paired DEG test.

\b
See also:
  Related:
    scrnah5adrawcounts2pseudobulk
  Depends:
    Python/Scanpy

\b
Date: 2023/03/14
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/python/{scriptname}.py'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"key={list(key)}",
		f"idheader='{idheader}'",
		f"separator='{separator}'",
		f"filterlabel='{filterlabel if filterlabel is not None else ''}'",
		f"mincell={mincell}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
