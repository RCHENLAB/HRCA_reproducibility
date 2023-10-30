#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exeR.exeR import exeR
from pathlib import Path
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-s', '--sampleid', type=click.STRING, default='sampleid', show_default=True, help='Label of sample ID.')
@click.option('-u', '--studyid', type=click.STRING, default='dataset', show_default=True, help='Label of study ID.')
@click.option('-c', '--celltype', type=click.STRING, default='majorclass', show_default=True, help='Label of cell type.')
@click.option('-W', '--width', type=click.FLOAT, default=6, show_default=True, help='Width.')
@click.option('-H', '--height', type=click.FLOAT, default=6, show_default=True, help='Height.')
@click.option('-m', '--margin', type=click.FLOAT, default=5, show_default=True, help='Margin for labels.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, sampleid, studyid, celltype, width, height, margin, infile):
	"""
To generate cell type similarity of study (or dataset).

INFILE is a .h5ad file.

\b
Example:
  f=$(parentsearch.sh -d test_data/scrnah5adhvg2metaneighborus hsfull152.h5ad)
  bname=$(basename "$f" .h5ad)
  outdir=$(mktemp -d -u)
  scrnah5adhvg2metaneighborus -d "$outdir" -b "$bname" -s sampleid -u dataset -c majorclass -W 6 -H 6 -m 5 -- "$f"

\b
See also:
  Depend:
    R/MetaNeighbor

\b
Date: 2023/05/08
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"sampleid='{sampleid}'",
		f"studyid='{studyid}'",
		f"celltype='{celltype}'",
		f"width={width}",
		f"height={height}",
		f"margin={margin}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	# os.chdir(outdir)
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
