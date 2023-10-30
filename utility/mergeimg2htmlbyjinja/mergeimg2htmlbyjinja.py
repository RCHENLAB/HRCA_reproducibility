#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import os
from pathlib import Path
import mimetypes
import base64
from jinja2 import Environment, FileSystemLoader
import click

# From snakemake/report/__init__.py
def data_uri(data, filename, encoding="utf8", mime="text/plain"):
	"""Craft a base64 data URI from file with proper encoding and mimetype."""
	data=base64.b64encode(data)
	uri="data:{mime};charset={charset};filename={filename};base64,{data}" "".format(
		filename=filename, mime=mime, charset=encoding, data=data.decode("utf-8")
		)
	return uri

def mime_from_file(file):
	mime, encoding=mimetypes.guess_type(file)
	if mime is None:
		mime="text/plain"
		print(f"Could not detect mimetype for {file}, assuming text/plain.", file=sys.stderr)
	return mime, encoding

def data_uri_from_file(file, defaultenc="utf8"):
	"""Craft a base64 data URI from file with proper encoding and mimetype."""
	if isinstance(file, Path):
		file=str(file)
	mime, encoding=mime_from_file(file)
	if encoding is None:
		encoding=defaultenc
	with open(file, "rb") as f:
		return data_uri(f.read(), os.path.basename(file), encoding, mime)

# From snakemake/report/data/common.py
def get_resource_as_string(path):
	return open(Path(__file__).parent / "template" / path).read()

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-n', '--headname', multiple=True, type=click.STRING, help='Header names for infile. Default: basenames of infiles.')
@click.option('-t', '--typename', multiple=True, type=click.STRING, help='Image types, available png or pdf for each file. Default: suffix of infiles.')
@click.argument('infile', type=click.Path(exists=True), nargs=-1)
def main(outdir, bname, headname, typename, infile):
	"""
Merge PNG/PDF files to a .html file using Jinja2 template.

\b
Example:
  indir=$(parentsearch.sh -d test_data mergepng2htmlbyjinja)
  files=(
  "$indir"/umapchen82_hackney_0.01_0.01_0.01_umap_leiden_0_wolabel.png
  "$indir"/umapchen82_hackney_0.05_0.01_0.01_umap_leiden_0_wolabel.png
  "$indir"/umapchen82_hackney_0.1_0.01_0.01_umap_leiden_0_wolabel.png
  "$indir"/umapchen82_hackney_0.5_0.01_0.01_umap_leiden_0_wolabel.png
  "$indir"/umapchen82_hackney_1_0.01_0.01_umap_leiden_0_wolabel.png
  "$indir"/umapchen82_hackney_1.5_0.01_0.01_umap_leiden_0_wolabel.png
  "$indir"/umapchen82_hackney_2_0.01_0.01_umap_leiden_0_wolabel.png
  )
  outdir=$(mktemp -d -u)
  bname=umap_level0
  mergeimg2htmlbyjinja -d "$outdir" -b "$bname" -- "${files[@]}"

\b
  # -n|--headname
  headname=($(basename.sh -s .png -- "${files[@]}" | sed -e "s/umapchen82_hackney_//" -e "s/_umap_leiden_0_wolabel//"))
  mergeimg2htmlbyjinja -d "$outdir" -b "$bname" $(basharr2cmdopts.sh -o -n -- "${headname[@]}") -- "${files[@]}"

\b
  # PDF files
  mergeimg2htmlbyjinja -d "$outdir" -b "$bname" -- *.pdf

\b
Note:
  1. File extensions are used to determine the HTML tags: .pdf for <iframe> and .png for <img>.
  2. Enable empty input for a stable execution.
  3. oldname: mergepng2htmlbyjinja for PNG files only. @2022/12/20.

\b
Date: 2023/01/06
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	env=Environment(loader=FileSystemLoader(Path(__file__).parent / "template"))
	env.filters["get_resource_as_string"]=get_resource_as_string
	template=env.get_template("index.html.jinja2")

	## data_uri of files
	if len(headname)!=len(infile):
		headname=[Path(file).stem for file in infile]
	if len(typename)!=len(infile):
		typename=[Path(file).suffix for file in infile]

	_results=[(name, type, data_uri_from_file(file)) for name, type, file in zip(headname, typename, infile)]

	## Render report file from a Jinja template
	Path(outdir).mkdir(parents=True, exist_ok=True)
	with open(f'{outdir}/{bname}.html', mode='w', encoding='utf-8') as f:
		f.write(
			template.render(
				results=_results,
				)
			)
	return 0

if __name__ == "__main__":
	sys.exit(main())
