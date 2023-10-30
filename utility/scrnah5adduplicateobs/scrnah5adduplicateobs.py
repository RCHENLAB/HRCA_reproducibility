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
@click.option('-s', '--src', multiple=True, type=click.STRING, required=True, help='Source metadata name.')
@click.option('-t', '--dest', multiple=True, type=click.STRING, required=True, help='Destination metadata name.')
@click.argument('infile', type=click.Path(exists=True))
def main(outdir, bname, condaenv, src, dest, infile):
	"""
Duplicate existing metadata names to new names.

INFILE is a .h5ad file.

\b
Example:
  slurmtaco.sh -t 2 -m 20G -- scrnah5adduplicateobs -d "$outdir" -b "$bname" -s leiden_1 -t saturn_label -- "$f"

\b
See also:
  Related:
    scrnah5adrenameobs.sh

\b
Date: 2023/02/07
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/python/{scriptname}.py'

	if len(src)!=len(dest):
		click.echo(f'{scriptname}: -s|--src and -t|--dest must have the same numbers of elements.', file=sys.stderr)
		sys.exit(1)

	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"src={src}",
		f"dest={dest}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
