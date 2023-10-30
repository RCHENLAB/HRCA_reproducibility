#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
import shutil
from exeR.exeR import exeR
from pathlib import Path
import click
from multiprocessing import Pool
from itertools import repeat

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-t', '--numthreads', type=click.INT, default=os.cpu_count(), show_default=True, help='Number of threads.')
@click.option('-b', '--bname', multiple=True, type=click.STRING, help='Bnames for infile. Default: basenames of infiles.')
@click.option('-W', '--width', type=click.FLOAT, default=1000, show_default=True, help='Width in pixel.')
@click.option('-H', '--height', type=click.FLOAT, default=1000, show_default=True, help='Height in pixel.')
@click.option('-z', '--zoom', type=click.INT, default=5, show_default=True, help='Zoom factor.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True), nargs=-1)
def main(outdir, condaenv, numthreads, bname, width, height, zoom, infile):
	"""
To screenshot HTML to PNG.

INFILE are HTML file[s]. Javascript is supported in the HTML files.

\b
Example:
  infile=$(parentsearch.sh -d test_data/html2png RGC_macaque_human_macaque_GradientBoostingClassifier_0.html)
  outdir=$(mktemp -d -u)
  html2png -d "$outdir" -- "$infile"

\b
  # An example with extended globbing in Bash
  shopt -s extglob
  html2png -W 500 -H 600 ../BC_*KNeighborsClassifier_@(0.9|0.99)_!(*max).html

\b
Note:
  1. The zoom factor (`-z|--zoom`) will result in a factor as many pixels vertically and horizontally.

\b
See also:
  Depends:
    R/webshot2::webshot()
    R/magick::image_trim()

\b
Date: 2023/02/17
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if len(bname)!=len(infile):
		# bname=[Path(file).name.removesuffix(''.join(Path(file).suffixes)) for file in infile]
		bname=[Path(file).stem for file in infile]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	with Pool(processes=numthreads) as p:
		p.starmap(html2png, zip(bname, infile, repeat(width), repeat(height), repeat(zoom), repeat(condaenv)))
	return 0

def html2png(*args):
	bname, infile, width, height, zoom, condaenv=args
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	exprs=[
		f"bname='{bname}'",
		f"width={width}",
		f"height={height}",
		f"zoom={zoom}",
		f"infile='{infile}'",
	]
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
