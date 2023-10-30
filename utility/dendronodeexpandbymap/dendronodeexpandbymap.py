#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from pathlib import Path
from multiprocess import Pool
import re
import pandas as pd
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-t', '--numthreads', type=click.INT, default=os.cpu_count(), show_default=True, help='Number of threads.')
@click.option('-b', '--bname', multiple=True, type=click.STRING, help='Bnames for infile. Default: basenames of infiles.')
@click.option('-m', '--nodemap', type=click.Path(exists=True), required=True, help='Node mapping.')
@click.option('-l', '--edgelen', type=click.FLOAT, default=0.1, show_default=True, help='Edge length.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True), nargs=-1)
def main(outdir, condaenv, numthreads, bname, nodemap, edgelen, infile):
	"""
To expand nodes from newick tree[s].

INFILE is newick tree file[s].

\b
Example:
  outdir=$(mktemp -d -u)
  mmfile=$outdir/input/mmfile.tree
  hsfile=$outdir/input/hsfile.tree
  mafile=$outdir/input/mafile.tree
  nodemap=$outdir/input/nodemap.txt
  tofile.sh -o "$mmfile" <<EOF
  (RBC:2.624137,((((BC5B:2.624137,BC5C:2.624137):5.094375,(BC5A:2.624137,BC5D:2.624137):3.702513):3.183063,((BC8:2.624137,BC9:2.624137):3.395735,(BC6:2.624137,BC7:2.624137):0.6367871):2.016773):0.6408071,(BC3A:2.624137,((BC3B:2.624137,BC4:2.624137):5.089302,(BC2:2.624137,(BC1A:2.624137,BC1B:2.624137):6.916711):3.770641):0.9821529):0.1895527):2.792187);
  EOF
  tofile.sh -o "$mafile" <<EOF
  (RB:1.67966,(DB6:1.67966,(DB3a:1.67966,((DB1:1.67966,OFFx:1.67966):1.804481,((FMB:1.67966,(DB2:1.67966,DB3b:1.67966):2.677551):2.488999,((DB4:1.67966,DB5*:1.67966):2.72492,(BB/GB*:1.67966,IMB:1.67966):2.707165):0.607485):0.9986878):1.131455):0.3387803):0.6710486);
  EOF
  tofile.sh -o "$hsfile" <<EOF
  (RBC:2.872565,(DB3a:2.872565,((DB6:2.872565,(DB1:2.872565,(FMB:2.872565,OFFx:2.872565):1.393701):3.746643):1.462098,((DB2:2.872565,DB3b:2.872565):9.192604,((BB/GB:2.872565,IMB:2.872565):3.801548,(DB5:2.872565,(DB4a:2.872565,DB4b:2.872565):5.653498):2.752362):1.397211):0.2335877):1.058024):3.402491);
  EOF
  tofile.sh -o "$nodemap" <<EOF
  BB/GB	hBB/GB
  BB/GB	aBB/GB
  BB/GB	mBC8/BC9
  RBC	hRB
  RBC	aRB
  RBC	mRBC
  EOF
  dendronodeexpandbymap -d "$outdir/output" -m "$nodemap" -- "$mmfile" "$hsfile" "$mafile"

\b
Note:
  1. A node map consists of two columns for key-value pairs. If a key has only one value, the 'edgelen' (specified by '-l|--edgelen') will not be used, and only the node label will be renamed. However, if a key has multiple values, one node will be expanded into several nodes to accommodate these value pairs, with an additional 'edgelen'.

\b
See also:
  Related:
    dendronewick2pdf

\b
Date: 2023/03/31
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if len(bname)!=len(infile):
		bname=[Path(file).stem for file in infile]

	# compile the pattern
	def col2join(x):
		if x.shape[0]>1: # expand node into multiple nodes with an edgelen.
			xstr=(x[1].astype(str)+f":{edgelen}").tolist()
			xstr=f"({','.join(xstr)})"
		else:
			xstr=x.iloc[0, 1] # rename the node label
		return xstr

	nodemap=pd.read_csv(nodemap, header=None, sep='\t')
	nodemap=nodemap.groupby(0).apply(col2join).to_dict()
	pattern=re.compile("|".join([rf"\b{k}\b" for k in nodemap.keys()]), flags=re.MULTILINE)

	# substitute
	def cmd(*args):
		name, file=args
		outfile=f"{outdir}/{name}.tree"
		with open(file, 'r') as reader:
			result=pattern.sub(lambda match: nodemap[match.group(0)], reader.read())
		with open(outfile, 'w') as writer:
			writer.write(result)

	Path(outdir).mkdir(parents=True, exist_ok=True)
	with Pool(processes=numthreads) as p:
		p.starmap(cmd, zip(bname, infile))
	return 0

if __name__ == "__main__":
	sys.exit(main())
