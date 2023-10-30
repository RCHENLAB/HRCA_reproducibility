#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exeR.exeR import exeR
from pathlib import Path
from multiprocess import Pool
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-t', '--numthreads', type=click.INT, default=os.cpu_count(), show_default=True, help='Number of threads.')
@click.option('-b', '--bname', multiple=True, type=click.STRING, help='Bnames for infile. Default: basenames of infiles.')
@click.option('-H', '--height', type=click.FLOAT, default=5, show_default=True, help='Height.')
@click.option('-W', '--width', type=click.FLOAT, default=5, show_default=True, help='Width.')
@click.option('-c', '--colormap', type=click.Path(exists=True, resolve_path=True), help='A 2-column headed color map file.')
@click.option('-r', '--direction', type=click.Choice(['rightwards', 'leftwards', 'upwards', 'downwards']), is_flag=False, flag_value='downwards', default='rightwards', show_default=True, help='Direction.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True), nargs=-1)
def main(outdir, condaenv, numthreads, bname, height, width, colormap, direction, infile):
	"""
To generate a PDF file from a newick file.

INFILE is a newick tree file.

\b
Example:
  mmfile=$(mktemp -u --suffix=.tree)
  hsfile=$(mktemp -u --suffix=.tree)
  mafile=$(mktemp -u --suffix=.tree)
  tofile.sh -o "$mmfile" <<EOF
  (RBC:2.624137,((((BC5B:2.624137,BC5C:2.624137):5.094375,(BC5A:2.624137,BC5D:2.624137):3.702513):3.183063,((BC8:2.624137,BC9:2.624137):3.395735,(BC6:2.624137,BC7:2.624137):0.6367871):2.016773):0.6408071,(BC3A:2.624137,((BC3B:2.624137,BC4:2.624137):5.089302,(BC2:2.624137,(BC1A:2.624137,BC1B:2.624137):6.916711):3.770641):0.9821529):0.1895527):2.792187);
  EOF
  tofile.sh -o "$mafile" <<EOF
  (RB:1.67966,(DB6:1.67966,(DB3a:1.67966,((DB1:1.67966,OFFx:1.67966):1.804481,((FMB:1.67966,(DB2:1.67966,DB3b:1.67966):2.677551):2.488999,((DB4:1.67966,'DB5*':1.67966):2.72492,('BB/GB*':1.67966,IMB:1.67966):2.707165):0.607485):0.9986878):1.131455):0.3387803):0.6710486);
  EOF
  tofile.sh -o "$hsfile" <<EOF
  (RBC:2.872565,(DB3a:2.872565,((DB6:2.872565,(DB1:2.872565,(FMB:2.872565,OFFx:2.872565):1.393701):3.746643):1.462098,((DB2:2.872565,DB3b:2.872565):9.192604,(('BB/GB':2.872565,IMB:2.872565):3.801548,(DB5:2.872565,(DB4a:2.872565,DB4b:2.872565):5.653498):2.752362):1.397211):0.2335877):1.058024):3.402491);
  EOF
  outdir=$(mktemp -d -u)
  dendronewick2pdf -d "$outdir" -H 5 -W 4 -- "$mmfile" "$hsfile" "$mafile"
  dendronewick2pdf -d "$outdir" -H 6 -W 5 -c "$colormap" -- "$treefile"

\b
Note:
  1. The color map is headed, e.g., label<TAB>species, and the second column will be the legend. E.g.,

\b
  ```
  $ head.sh colormap.txt
  ==> colormap.txt <==
  label	species
  hBB/GB_1	human
  hBB/GB_2	human
  ```

\b
See also:
  Depends:
    R/phylogram

\b
Date: 2023/03/23
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if len(bname)!=len(infile):
		bname=[Path(file).stem for file in infile]
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	Path(outdir).mkdir(parents=True, exist_ok=True)
	# os.chdir(outdir)

	def newick2pdf(*args):
		name, file=args
		exprs=[
			f"outdir='{outdir}'",
			f"bname='{name}'",
			f"height={height}",
			f"width={width}",
			f"colormap='{colormap if colormap is not None else ''}'",
			f"direction='{direction}'",
			f"infile='{file}'",
			]
		exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

	with Pool(processes=numthreads) as p:
		p.starmap(newick2pdf, zip(bname, infile))
	return 0

if __name__ == "__main__":
	sys.exit(main())
