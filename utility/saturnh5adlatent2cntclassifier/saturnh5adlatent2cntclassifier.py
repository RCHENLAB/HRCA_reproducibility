#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exePython.exePython import exePython
from pathlib import Path
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)
def main():
	"""
Develop conventional classifiers using SATURN latent representations.

\b
Example:
  saturnh5adlatent2cntclassifier -h
  saturnh5adlatent2cntclassifier list
  saturnh5adlatent2cntclassifier cv -h
  saturnh5adlatent2cntclassifier train -h
  saturnh5adlatent2cntclassifier trainby -h
  saturnh5adlatent2cntclassifier predict -h
  saturnh5adlatent2cntclassifier predictproba -h

\b
Date: 2023/02/07
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	pass

# List models
@main.command(context_settings=CONTEXT_SETTINGS)
def list():
	"""
List implemented models.

\b
Note:
  1. Skip `GradientBoostingClassifier` and `NeighborhoodComponentsAnalysis` due to too slow for training. @2/16/2023

\b
Date: 2023/02/16
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	models=[
		"MultinomialNB",
		"SVC",
		"KNeighborsClassifier",
		"LogisticRegression",
		"RandomForestClassifier",
		"SVClinear",
		"SGDClassifier",
		"NearestCentroid",
		"RadiusNeighborsClassifier",
		# "GradientBoostingClassifier",
		# "NeighborhoodComponentsAnalysis",
		]
	click.echo('\n'.join(models))
	return 0

# Cross-validation
@main.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-m', '--model', type=click.STRING, required=True, help='Conventional multi-class classification model.')
@click.option('-p', '--params', type=click.STRING, help='Keyword arguments in a dict string.')
@click.option('-c', '--classlabel', type=click.STRING, default='celltype', show_default=True, help='Classification label.')
@click.option('-s', '--seed', type=click.INT, default=12345, show_default=True, help='Random seed.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def cv(outdir, bname, condaenv, model, params, classlabel, seed, infile):
	"""
10-fold cross validate one classifier by model.

INFILE is a .h5ad file.

\b
Example:
  indir=$(mrrdir.sh ../../scrnah5adsplitby)
  outdir=$(mrrdir.sh)
  function cmd {
  local f=$1
  local model=$2
  local bname=$(basename "$f" .h5ad)_$model
  if fileexists.sh "$f"
  then
  	slurmtaco.sh -t 2 -m 20G -- saturnh5adlatent2cntclassifier cv -d "$outdir" -b "$bname" -m "$model" -c labels2 -e saturn -- "$f"
  	slurmtaco.sh -t 2 -m 20G -- saturnh5adlatent2cntclassifier cv -d "$outdir" -b "$bname" -m "$model" -p "{}" -c labels2 -e saturn -- "$f"
  fi
  }
  source env_parallel.bash
  env_parallel cmd ::: "$indir"/test256_data_BC_mm_concat47_inner_org_saturn_seed_12345_mouse.h5ad ::: $(saturnh5adlatent2cntclassifier list)

\b
See also:
  Upstream:
    saturntraincrossspecies
    scrnah5adsplitby
  Related:
    merfishh5ad2cntclassifier10cv.sh
    merfishh5ad2cntclassifiertrain.sh
    merfishh5ad2cntclassifierpredict.sh

\b
Date: 2023/02/07
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f"{absdir}/python/cv.py"
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"model='{model}'",
		f"params={params if params is not None else {}}",
		f"classlabel='{classlabel}'",
		f"seed={seed}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

# Train a classifier
@main.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-m', '--model', type=click.STRING, required=True, help='Conventional multi-class classification model.')
@click.option('-c', '--classlabel', type=click.STRING, default='celltype', show_default=True, help='Classification label.')
@click.option('-s', '--seed', type=click.INT, default=12345, show_default=True, help='Random seed.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def train(outdir, bname, condaenv, model, classlabel, seed, infile):
	"""
Train a classifier using the full data.

INFILE is a .h5ad file.

\b
Example:
  indir=$(mrrdir.sh ../../scrnah5adsplitby)
  outdir=$(mrrdir.sh)
  function cmd {
  local f=$1
  local model=$2
  local bname=$(basename "$f" .h5ad)_$model
  if fileexists.sh "$f"
  then
  	slurmtaco.sh -t 2 -m 20G -- saturnh5adlatent2cntclassifier train -d "$outdir" -b "$bname" -m "$model" -c labels2 -e saturn -- "$f"
  fi
  }
  source env_parallel.bash
  env_parallel cmd ::: "$indir"/test256_data_BC_mm_concat47_inner_org_saturn_seed_12345_mouse.h5ad ::: $(saturnh5adlatent2cntclassifier list)

\b
See also:
  Upstream:
    saturntraincrossspecies
    scrnah5adsplitby
  Related:
    merfishh5ad2cntclassifier10cv.sh
    merfishh5ad2cntclassifiertrain.sh
    merfishh5ad2cntclassifierpredict.sh

\b
Date: 2023/02/07
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f"{absdir}/python/train.py"
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"model='{model}'",
		f"classlabel='{classlabel}'",
		f"seed={seed}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

# Train a classifier by customized parameters
@main.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-m', '--model', type=click.STRING, required=True, help='Conventional multi-class classification model.')
@click.option('-p', '--params', type=click.STRING, help='Keyword arguments in a dict string.')
@click.option('-c', '--classlabel', type=click.STRING, default='celltype', show_default=True, help='Classification label.')
@click.option('-s', '--seed', type=click.INT, default=12345, show_default=True, help='Random seed.')
@click.option('-k', '--skipscale', is_flag=True, help='Skip StandardScaler().')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def trainby(outdir, bname, condaenv, model, params, classlabel, seed, skipscale, infile):
	"""
Train a classifier by customized parameters.

INFILE is a .h5ad file.

\b
Example:
  f=$(parentsearch.sh -d test_data/saturnh5adlatent2cntclassifier BC_mouse_macaque_human.h5ad)
  outdir=$(mktemp -d -u)
  bname=$(basename "$f" .h5ad)
  saturnh5adlatent2cntclassifier trainby -c labels -s 12345 -e scanpyanndata -d "$outdir" -b "$bname" -m "KNeighborsClassifier" -p "{'n_neighbors': 100}" -- "$f"

\b
Note:
  1. `-p|--params` is a good feature to support variable arguments for the function. Empty directionary is allowed.

\b
See also:
  Upstream:
    saturntraincrossspecies
    scrnah5adsplitby
  Related:
    merfishh5ad2cntclassifier10cv.sh
    merfishh5ad2cntclassifiertrain.sh
    merfishh5ad2cntclassifierpredict.sh

\b
Date: 2023/02/07
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f"{absdir}/python/trainby.py"
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"model='{model}'",
		f"params={params if params is not None else {}}",
		f"classlabel='{classlabel}'",
		f"seed={seed}",
		f"skipscale={skipscale}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

# Predict labels by a trained classifier
@main.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-m', '--model', type=click.Path(exists=True, resolve_path=True), required=True, help='Conventional multi-class classification model.')
@click.option('-l', '--label', type=click.STRING, help='A column label to evaluate prediction when available.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def predict(outdir, bname, condaenv, model, label, infile):
	"""
Predict labels by a trained classifier.

INFILE is a .h5ad file.

\b
Example:
  indir=$(mrrdir.sh ../../scrnah5adsplitby)
  modeldir=$(mrrdir.sh ../train)
  outdir=$(mrrdir.sh)
  function cmd {
  local f=$1
  local modelname=$2
  local modelfile=$modeldir/test256_data_BC_mm_concat47_inner_org_saturn_seed_12345_mouse_${modelname}_classifier.pbz2
  local bname=$(basename "$f" .h5ad)_$modelname
  if fileexists.sh "$f"
  then
  	slurmtaco.sh -t 2 -m 20G -- saturnh5adlatent2cntclassifier predict -d "$outdir" -b "$bname" -m "$modelfile" -e saturn -- "$f"
  	slurmtaco.sh -t 2 -m 20G -- saturnh5adlatent2cntclassifier predict -d "$outdir" -b "$bname" -m "$modelfile" -e saturn -l labels2 -- "$f"
  fi
  }
  source env_parallel.bash
  env_parallel cmd ::: "$indir"/test256_data_BC_mm_concat47_inner_org_saturn_seed_12345_human.h5ad ::: $(saturnh5adlatent2cntclassifier list)

\b
Note:
  1. `-l|--label` might not match the predicted label in the trained classifier. So, metrics might not be valid, except ARI.

\b
See also:
  Upstream:
    saturnh5adlatent2cntclassifier list
    saturnh5adlatent2cntclassifier train
  Related:
    merfishh5ad2cntclassifier10cv.sh
    merfishh5ad2cntclassifiertrain.sh
    merfishh5ad2cntclassifierpredict.sh

\b
Date: 2023/02/07
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f"{absdir}/python/predict.py"
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"model='{model}'",
		f"label='{label if label is not None else ''}'",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

# Predict labels by a trained classifier with a threshold
@main.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-m', '--model', type=click.Path(exists=True, resolve_path=True), required=True, help='Conventional multi-class classification model.')
@click.option('-l', '--label', type=click.STRING, help='A column label to evaluate prediction when available.')
@click.option('-r', '--threshold', type=click.FLOAT, default=0.9, show_default=True, help='A threshold of classification probabilities.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def predictproba(outdir, bname, condaenv, model, label, threshold, infile):
	"""
Predict labels by a threshold of classification probabilities.

INFILE is a .h5ad file.

\b
Example:
  indir=$(mrrdir.sh ../../scrnah5adsplitby)
  modeldir=$(mrrdir.sh ../train)
  outdir=$(mrrdir.sh)
  function cmd {
  local f=$1
  local modelname=$2
  local modelfile=$modeldir/test256_data_BC_mm_concat47_inner_org_saturn_seed_12345_mouse_${modelname}_classifier.pbz2
  local bname=$(basename "$f" .h5ad)_$modelname
  if fileexists.sh "$f"
  then
  	slurmtaco.sh -t 2 -m 20G -- saturnh5adlatent2cntclassifier predictproba -d "$outdir" -b "$bname" -m "$modelfile" -e saturn -- "$f"
  	slurmtaco.sh -t 2 -m 20G -- saturnh5adlatent2cntclassifier predictproba -d "$outdir" -b "$bname" -m "$modelfile" -e saturn -l labels2 -- "$f"
  fi
  }
  source env_parallel.bash
  env_parallel cmd ::: "$indir"/test256_data_BC_mm_concat47_inner_org_saturn_seed_12345_human.h5ad ::: $(saturnh5adlatent2cntclassifier list)

\b
Note:
  1. `-l|--label` might not match the predicted label in the trained classifier. So, metrics might not be valid, except ARI.

\b
See also:
  Upstream:
    saturnh5adlatent2cntclassifier list
    saturnh5adlatent2cntclassifier train
  Related:
    saturnh5adlatent2cntclassifier predict
    merfishh5ad2cntclassifier10cv.sh
    merfishh5ad2cntclassifiertrain.sh
    merfishh5ad2cntclassifierpredict.sh

\b
Date: 2023/02/28
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f"{absdir}/python/predictproba.py"
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"model='{model}'",
		f"label='{label if label is not None else ''}'",
		f"threshold={threshold}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

# Main interface
if __name__ == "__main__":
	sys.exit(main())
