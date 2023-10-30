#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
from pathlib import Path
import pandas as pd
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-x', '--xlabel', type=click.STRING, multiple=True, default=['saturn_label'], show_default=True, help='X-Labels.')
@click.option('-y', '--ylabel', type=click.STRING, multiple=True, default=['predict'], show_default=True, help='X-Labels.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, xlabel, ylabel, infile):
	"""
Tabulate occurences of two metadata from a TSV file.

\b
Example:
  f=$(parentsearch.sh -d test_data/saturncntclassifierpred2tab BC_human_mouse_LogisticRegression.txt.gz) 
  bname=$(basename "$f" .txt.gz)
  outdir=$(mktemp -d -u)
  head.sh "$f"
  xlabels=(
  group1
  group2
  group3
  cluster1
  cluster2
  )
  ylabels=(
  predict
  )
  saturncntclassifierpred2tab -d "$outdir" -b "$bname" $(basharr2cmdopts.sh -o -x -- "${xlabels[@]}") $(basharr2cmdopts.sh -o -y -- "${ylabels[@]}") -- "$f"

\b
Note:
  1. A label with < 2 categories will be skipped.

\b
See also:
  Upstream:
    tsvfileaddfromcolumn
  Downstream:
    rtable2sankeydiagram

\b
Date: 2023/02/23
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	input=pd.read_csv(infile, sep='\t', header=0) if infile.endswith(('.txt', '.txt.gz', '.tsv', 'tsv.gz')) else pd.read_csv(infile, header=0)
	Path(outdir).mkdir(parents=True, exist_ok=True)
	for xl, yl in ((x, y) for x in xlabel for y in ylabel):
		if all(h in input.columns for h in [xl, yl]):
			if len(input[xl].unique())>1 and len(input[yl].unique())>1:
				tab=pd.crosstab(input[xl], input[yl])
				tab.to_csv(f"{outdir}/{bname}_{xl}_{yl}.txt.gz", sep='\t')
			else:
				click.echo(f"Warning: {xl} or {yl} have less than 2 categories.")
		else:
			click.echo(f"Warning: {xl} or {yl} is not in the header.")
	return 0

if __name__ == "__main__":
	sys.exit(main())
