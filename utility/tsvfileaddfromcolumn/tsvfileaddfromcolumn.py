#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import pandas as pd
from pathlib import Path
import click
CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=False, help='Bname for an outfile. Print to stdout by default.')
@click.option('-s', '--srcfile', type=click.Path(exists=True), required=True, help='Source file.')
@click.option('-k', '--key', type=click.STRING, multiple=True, help='Keys to extract columns. Full columns by default.')
@click.argument('infile', type=click.Path(exists=True))
def main(outdir, bname, srcfile, key, infile):
	"""
Append columns from a file.

\b
Example:
  indir=$(mktemp -d -u)
  file1=$indir/file1.txt.gz
  file2=$indir/file2.csv
  tofile.sh -o "$file1" <<EOF
  c1	c2
  a	1
  b	2
  EOF
  tofile.sh -o "$file2" <<EOF
  c1,c3
  c,x
  d,y
  EOF
  outdir=$indir
  bname=output
  tsvfileaddfromcolumn -s "$file2" -- "$file1"
  tsvfileaddfromcolumn -s "$file2" -k c3 -- "$file1"
  tsvfileaddfromcolumn -d "$outdir" -b "$bname" -s "$file2" -k c1 -k c3 -- "$file1"
  # bug fix notin keys
  tsvfileaddfromcolumn -s "$file2" -k notin -k c2 -k c1 -- "$file1"

\b
Note:
  1. File type will be inferred from the extension.

\b
Date: 2023/02/21
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	input, src=[
		pd.read_csv(file, sep='\t', header=0)
		if file.endswith(('.txt', '.txt.gz', '.tsv', 'tsv.gz'))
		else pd.read_csv(file, header=0)
		for file in [infile, srcfile]
		]
	if key:
		src=src.loc[:, src.columns.isin(key)].copy()

	input=pd.concat([input, src], axis=1)
	if bname:
		Path(outdir).mkdir(parents=True, exist_ok=True)
		input.to_csv(f'{outdir}/{bname}.txt.gz', sep='\t', index=False)
	else:
		input.to_csv(sys.stdout, sep='\t', index=False)
	return 0

if __name__ == "__main__":
	sys.exit(main())
