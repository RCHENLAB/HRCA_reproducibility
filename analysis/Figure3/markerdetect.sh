#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug

function cmd {
local f=$1
local model=$2
local ncombn=$3
local ntop=$4
local bname=$(basename "$f" .h5ad)_ntop${ntop}_${model}_${ncombn}

scrnah5ad2markerbinaryaurochvg -d "$outdir" -b "$bname" -l cluster2 -m "$model" -k sampleid -H 10000 -n "$ntop" -c "$ncombn" -t 8 -N -- "$f"
}

source env_parallel.bash
env_parallel cmd ::: snRNA_AC_rawcount.h5ad ::: LogisticRegression RandomForestClassifier SVC XGBClassifier ::: 1 2 3 ::: 10 20 50
env_parallel cmd ::: snRNA_RGC_rawcount.h5ad ::: LogisticRegression RandomForestClassifier SVC XGBClassifier ::: 1 2 3 ::: 10 20 50
