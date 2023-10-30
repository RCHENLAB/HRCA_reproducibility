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
@click.option('-t', '--iteration', type=click.INT, default=32, show_default=True, help='Iterations to optimize diagram layout.')
@click.option('-s', '--saveimage', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Save the screenshot to a PNG file. (Chrome should be installed).')
@click.option('-W', '--width', type=click.FLOAT, default=1000, show_default=True, help='Width in pixel for `-s|--saveimage`.')
@click.option('-H', '--height', type=click.FLOAT, default=1000, show_default=True, help='Height in pixel.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, iteration, width, height, saveimage, infile):
	"""
To plot Sankey diagram from a contingency table.

INFILE is a .txt.gz file.

\b
Example:
  infile=$(mktemp -u --suffix=.txt.gz)
  # BC_mouse_macaque_human_mouse_LogisticRegression_tab.txt.gz
  tofile.sh -o "$infile" <<EOF
  labels2	BC1A	BC1B	BC2	BC3A	BC3B	BC4	BC5A	BC5C	BC5D	BC6	BC7	BC8	BC9	RBC
  HBC1	1999	0	0	0	0	0	0	0	0	0	0	0	0	0
  HBC10	0	0	0	0	0	0	0	0	0	1998	0	1	0	0
  HBC11	0	1999	0	0	0	0	0	0	0	0	0	0	0	0
  HBC12	0	0	0	1990	0	0	0	0	0	0	0	0	0	0
  HBC13	0	0	0	0	0	0	0	0	1993	0	0	0	0	0
  HBC14	0	0	0	0	0	1	0	0	0	0	0	37	1958	0
  HBC2	0	0	0	0	0	0	0	0	0	0	0	0	0	1995
  HBC3	0	0	0	0	0	0	0	0	0	0	1998	0	0	0
  HBC4	0	0	0	0	0	1998	0	0	0	0	0	0	0	0
  HBC5	0	0	0	0	0	0	0	2000	0	0	0	0	0	0
  HBC6	0	0	1998	0	0	0	0	0	0	0	0	0	0	0
  HBC7	0	0	0	0	0	0	1997	0	0	0	0	0	0	0
  HBC8	0	0	0	0	1998	0	0	0	0	0	0	0	0	0
  HBC9	0	0	0	0	0	0	0	0	0	0	0	1722	278	0
  EOF
  outdir=$(mktemp -d -u)
  bname=BC_hs_mm_LR
  rtable2sankeydiagram -d "$outdir" -b "$bname" -- "$infile"
  rtable2sankeydiagram -d "$outdir" -b "$bname" -t 0 -- "$infile"
  rtable2sankeydiagram -s -d "$outdir" -b "$bname" -t 0 -- "$infile"

\b
Note:
  1. `-t 0` will turn off the optimization of the diagram layout. The nodes are numerically sorted.
  2. Only limited color scheme is supported.
  https://github.com/christophergandrud/networkD3/issues/248#issuecomment-436968165

\b
  ```I.e.,
  d3.schemeCategory10
  d3.schemeCategory20
  d3.schemeCategory20b
  d3.schemeCategory20c
  ```

\b
  2.1 Other schemes from d3-scale-chromatic are not implemented in networkD3. So, use the default color scheme, i.e., `d3.schemeCategory20`.
  https://www.npmjs.com/package/d3-scale-chromatic

\b
  2.2 A customized color scaling

\b
  ```
  ColourScal ='d3.scaleOrdinal().range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'
  ```

\b
  3. Can run a local PC only. R/webshot2::webshot() needs Chrome installed on a local PC.

\b
  ```
  `google-chrome` and `chromium-browser` were not found. Try setting the CHROMOTE_CHROME environment variable or adding one of these executables to your PATH.
  Error in initialize(...) : Invalid path to Chrome
  Calls: source ... <Anonymous> -> initialize -> <Anonymous> -> initialize
  Execution halted
  ```

\b
See also:
  Upstream:
    saturnh5adlatent2cntclassifier
  Depends:
    R/networkD3
    R/webshot2

\b
Date: 2023/02/14
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"iteration={iteration}",
		f"saveimage={saveimage}",
		f"width={width}",
		f"height={height}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
