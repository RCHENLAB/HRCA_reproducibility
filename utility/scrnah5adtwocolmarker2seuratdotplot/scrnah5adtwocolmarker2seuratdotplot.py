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
@click.option('-W', '--width', type=click.FLOAT, default=8, show_default=True, help='Width.')
@click.option('-H', '--height', type=click.FLOAT, default=8, show_default=True, help='Height.')
@click.option('-x', '--xlab', type=click.STRING, help='X-axis label.')
@click.option('-y', '--ylab', type=click.STRING, help='Y-axis label.')
@click.option('-g', '--groupby', type=click.STRING, default='leiden', show_default=True, help='Metadata to group cells.')
@click.option('-m', '--marker', type=click.Path(exists=True), required=True, help='A two-column marker file.')
@click.option('-c', '--ordergroup', type=click.Path(exists=True), help='The order of the group labels.')
@click.option('-n', '--norm', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Normalize raw counts.')
@click.option('-f', '--flip', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Flip axis.')
@click.option('-F', '--bit64double', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Bit64 conversion to double for a big infile.')
@click.option('--scalemin', type=click.FLOAT, help='Lower limit for scaling.')
@click.option('--scalemax', type=click.FLOAT, help='Upper limit for scaling.')
@click.argument('infile', type=click.Path(exists=True))
def main(outdir, bname, condaenv, width, height, xlab, ylab, groupby, marker, ordergroup, norm, flip, bit64double, scalemin, scalemax, infile):
	"""
Generate a dotplot from a .h5ad file and an uheaded two-column marker file `group<TAB>marker`.

INFILE is a .h5ad file.

\b
Example:
  f=/tmp/scRNA_BC.h5ad
  bname=$(basename "$f" .h5ad)
  marker=$(mktemp -d -u --suffix=.txt.gz)
  tofile.sh -o "$marker" <<EOF
  HBC1	PTPRO
  HBC2	PRKCA
  HBC3	DLG2
  HBC4	DSCAM
  HBC5	PRKCA
  HBC6	FAM155A
  HBC7	KIRREL3
  HBC8	KIRREL3
  HBC9	MAP3K5
  HBC10	ROBO1
  HBC11	FAM155A
  HBC12	ERBB4
  HBC13	KIRREL3
  HBC14	AFF3
  EOF
  outdir=$(mktemp -d -u)
  scrnah5adtwocolmarker2seuratdotplot -d "$outdir" -b "$bname" -m "$marker" -g cluster2 -f -n -- "$f"

\b
See also:
  Depends:
    R/Seurat

\b
Date: 2023/05/18
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	exprs=[
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
		f"bit64double={bit64double}",
		f"scalemin={scalemin}" if scalemin else f"scalemin=NA",
		f"scalemax={scalemax}" if scalemax else f"scalemax=NA",
		f"infile='{infile}'",
		]
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
