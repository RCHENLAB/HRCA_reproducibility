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
@click.option('-W', '--width', type=click.FLOAT, default=8, show_default=True, help='Width.')
@click.option('-H', '--height', type=click.FLOAT, default=8, show_default=True, help='Height.')
@click.option('-x', '--xlab', type=click.STRING, help='X-axis label.')
@click.option('-y', '--ylab', type=click.STRING, help='Y-axis label.')
@click.option('-g', '--groupby', type=click.STRING, default='leiden', show_default=True, help='Metadata to group cells.')
@click.option('-m', '--marker', type=click.Path(exists=True), required=True, help='A file of marker genes for cell types.')
@click.option('-c', '--ordergroup', type=click.Path(exists=True), help='The order of the group labels.')
@click.option('-n', '--norm', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Normalize raw counts.')
@click.option('-f', '--flip', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Flip axis.')
@click.option('--scalemin', type=click.FLOAT, help='Lower limit for scaling.')
@click.option('--scalemax', type=click.FLOAT, help='Upper limit for scaling.')
@click.option('-F', '--bit64double', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Bit64 conversion to double for a big infile.')
@click.argument('infile', type=click.Path(exists=True))
def main(outdir, bname, width, height, xlab, ylab, groupby, marker, ordergroup, norm, flip, scalemin, scalemax, bit64double, infile):
	"""
Generate a dotplot using Seurat::DotPlot() from a .h5ad file.

INFILE is a Scanpy object in a .h5ad format.

\b
Example:
  f=/tmp/aaa/chen82_hackney.h5ad
  bname=$(basename "$f" .h5ad)
  marker=$(retinamarkermajortype.sh --hs --RGC -p)
  outdir=$(mktemp -d -u)
  scrnah5ad2seuratdotplot -n -d "$outdir" -b "$bname" -m "$marker" -- "$f"
  scrnah5ad2seuratdotplot -n -d "$outdir" -b "$bname" -m "$marker" -y leiden -g leiden -- "$f"
  # ungrouped features
  scrnah5ad2seuratdotplot -n -d "$outdir" -b "$bname" -m <(cat.sh -q "$marker" | cut -f 2) -y leiden -g leiden -- "$f"

\b
  # A second example
  f=/tmp/scRNA_BC.h5ad
  bname=$(basename "$f" .h5ad)
  marker=$(mktemp -u --suffix=.txt)
  retinamarkermajortype.sh --hs --BC  | cut -f 2 | tac | tofile.sh -o "$marker"
  outdir=$(mktemp -d -u)
  scrnah5ad2seuratdotplot -n -d "$outdir" -b "$bname" -m "$marker" -g group6 -c order.txt -f -H 4 -W 6 -- "$f"
  scrnah5ad2seuratdotplot -n -d "$outdir" -b "$bname" -m "$marker" -g group6 -c order.txt -f -H 4 -W 6 --scalemax 75 -- "$f"

\b
Note:
  1. `-m|--marker` is unheaded. E.g.,

\b
  ```
  $ retinamarkermajortype.sh --hs --RGC | head.sh
  pan-RGC	RBPMS
  pan-RGC	THY1
  MGC_ON	TPBG
  ```

\b
  2. `-n|--norm` is optional. Click@8.1.x is needed to assign `flag_value='T'` when `[F|T]` is not provided to `-n|--norm`.
  https://click.palletsprojects.com/en/8.1.x/options/#optional-value

\b
See also:
  Depends:
    R/Seurat

\b
Date: 2022/11/10
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	exprs=[
		f"f='{infile}'",
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"width={width}",
		f"height={height}",
		f"xlab='{xlab}'" if xlab else f"xlab=NULL",
		f"ylab='{ylab}'" if ylab else f"ylab=NULL",
		f"groupby='{groupby}'",
		f"markerfile='{marker}'",
		f"ordergroup='{ordergroup if ordergroup is not None else ''}'",
		f"norm={norm}",
		f"flip={flip}",
		f"scalemin={scalemin}" if scalemin else f"scalemin=NA",
		f"scalemax={scalemax}" if scalemax else f"scalemax=NA",
		f"bit64double={bit64double}",
		]
	return exeR.callback(exprs, script=script, condaenv=None, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
