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
@click.option('-l', '--label', type=click.STRING, default='cluster', show_default=True, help='Label of cell cluster.')
@click.option('-m', '--model', type=click.STRING, required=True, help='Conventional multi-class classification model.')
@click.option('-p', '--params', type=click.STRING, help='Keyword arguments in a dict string.')
@click.option('-s', '--seed', type=click.INT, default=12345, show_default=True, help='Random seed.')
@click.option('-k', '--batchkey', type=click.STRING, help='Batch key for HVG.')
@click.option('-H', '--nhvg', type=click.INT, default=10000, show_default=True, help='Number of HVGs for top genes.')
@click.option('-n', '--ngene', type=click.INT, default=20, show_default=True, help='Top genes per cluster.')
@click.option('-c', '--ncombn', type=click.INT, default=3, show_default=True, help='Number of genes for combinations.')
@click.option('-t', '--numthreads', type=click.INT, default=8, show_default=True, help='Number of threads.')
@click.option('-N', '--norm', is_flag=True, help='Normalize and log1p transform raw counts.')
@click.option('-f', '--topgenefile', type=click.Path(exists=True), help='To reuse top-genes file.')
@click.option('-g', '--group', type=click.STRING, multiple=True, help='Sub-group[s] for detection.')
@click.option('--log2foldchange', type=click.FLOAT, default=0, show_default=True, help='log2foldchange for significant records.')
@click.option('--qvalue', type=click.FLOAT, default=1, show_default=True, help='Q-value for significant records.')
@click.option('-v', '--verbose', is_flag=True, help='Verbose.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, label, model, params, seed, batchkey, nhvg, ngene, ncombn, numthreads, norm, topgenefile, group, log2foldchange, qvalue, verbose, infile):
	"""
To detect marker gene sets (2 or 3 genes) for cell clusters. Only HVGs are used to enumerate combinations of gene sets.

INFILE is a .h5ad file.

\b
Example:
  f=/storage/singlecell/jinli/wkfl/atlashumanprj/integration/snRNA/proc/BC/novel_marker/auroc/scrnah5adsubsetsamplingbykey/scrnah5adsubsetbyvaluecounts/snRNA_BC.h5ad
  outdir=$(mrrdir.sh)
  function cmd {
  local f=$1
  local bname=$(basename "$f" .h5ad)
  if fileexists.sh "$f"
  then
    slurmtaco.sh -t 2 -m 20G -n d03 -- scrnah5ad2markerbinaryaurochvg -d "$outdir" -b "$bname" -l cluster2 -m LogisticRegression -k sampleid -H 10000 -n 6 -c 3 -t 16 --norm -- "$f"
    slurmtacc.sh -t 1 -m 2G -p rtx -- "$(cat <<-EOF
      scrnah5ad2markerbinaryaurochvg -e xgboost -d "$outdir" -b "$bname" -l cluster2 -m XGBClassifier -p "{'tree_method': 'gpu_hist', 'gpu_id': 0}" -k sampleid -H 10000 -n 20 -c 3 -t 16 --norm -- "$f"
    EOF
    )"
  fi
  }
  source env_parallel.bash
  env_parallel cmd ::: "$f"

\b
Note:
  1. Methods
  1.1 Rationale is to predict cell clusters by marker genes directly.
  1.2 Normalization is needed, but scaling might not be.
  1.3 Balanced cells are necessary. See `scrnah5adsubsetsamplingbykey.sh`.
  1.4 AUROC metrics or F1-score
  1.5 How to resolve the sample batches?

\b
  2. Use case for XGBClassifier
  2.1 'gpu_id' is needed to specify the right device, otherwise, device 0 will be used.
  2.2 'n_estimators' is 100 by default.

\b
  3. Available model names:
    - XGBClassifier
    - LogisticRegression
    - RandomForestClassifier
    - SVC
    - SVClinear
    - GradientBoostingClassifier
    - KNeighborsClassifier
    - RadiusNeighborsClassifier
    - MultinomialNB
    - NeighborhoodComponentsAnalysis
    - SGDClassifier

\b
  4. This is to speed up by training binary classifiers only for a combination of gene sets from one class HVGs. A multiclass classification is unnecessary (see scrnah5ad2clsmarkerbyaurochvg).

\b
  5. `-N|--norm` is essential to achieve better classification.

\b
See also:
  Related:
    scrnah5ad2clsmarkerbyauroc
    scrnah5ad2clsmarkerbyaurochvg

\b
Date: 2023/05/03
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/python/{scriptname}.py'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"label='{label}'",
		f"model='{model}'",
		f"params={params if params else {}}",
		f"seed={seed}",
		f"batchkey='{batchkey}'" if batchkey else f"batchkey=None",
		f"nhvg={nhvg}",
		f"ngene={ngene}",
		f"ncombn={ncombn}",
		f"numthreads={numthreads}",
		f"norm={norm}",
		f"topgenefile='{topgenefile}'" if topgenefile else f"topgenefile=None",
		f"group={group}",
		f"log2foldchange={log2foldchange}",
		f"qvalue={qvalue}",
		f"verbose={verbose}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	# os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
