#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug

# 1. Subsampling
scrnah5adsubsetsamplingbykey.sh -d "$outdir" -b snRNA_AC -k cluster2 -n 500 -- snRNA_AC_rawcount.h5ad # AC
scrnah5adsubsetsamplingbykey.sh -d "$outdir" -b snRNA_RGC -k cluster2 -n 500 -- snRNA_RGC_rawcount.h5ad # RGC

# 2. marker detection for AC
function cmd {
local f=$1
local model=$2
local ncombn=$3
local ntop=$4
local bname=$(basename "$f" .h5ad)_ntop${ntop}_${model}_${ncombn}

scrnah5ad2markerbinaryaurochvg -d "$outdir" -b "$bname" -l cluster2 -m "$model" -k sampleid -H 10000 -n "$ntop" -c "$ncombn" -t 8 -N -- "$f"
}

source env_parallel.bash
env_parallel cmd ::: snRNA_AC.h5ad ::: LogisticRegression RandomForestClassifier SVC XGBClassifier ::: 1 2 3 ::: 10 20 50 # AC

# 3. marker detection for RGC
function cmd {
local f=$1
local ntop=$2
local model=$3
local ncombn=$4
local bname=$(basename "$f" .h5ad)_ntop${ntop}_${model}_${ncombn}
scrnah5ad2markerbinaryaurochvg -d "$outdir" -b "$bname" -l cluster2 -m "$model" -k sampleid -H 10000 -n "$ntop" -c "$ncombn" -t 8 -N --log2foldchange 1.0 --qvalue 0.05 -- "$f"
}

source env_parallel.bash
env_parallel --colsep='\t' cmd <<EOF
snRNA_RGC.h5ad	50	SVC	1
snRNA_RGC.h5ad	50	SVC	2
snRNA_RGC.h5ad	50	SVC	3
EOF


# 4. Dot plot for AC and RGC
scrnah5adtwocolmarker2seuratdotplot -d "$outdir" -b "$bname" -m snRNA_AC_ntop50_SVC_1_f1.txt.gz -g cluster2 -f -n -H 12 -W 15 -c marker_AC_order.txt -- snRNA_AC.h5ad
scrnah5adtwocolmarker2seuratdotplot -d "$outdir" -b "$bname" -m snRNA_RGC_marker.txt -g RGC_celltype -f -n -H 18 -W 10 -c marker_RGC_order.txt -- snRNA_RGC.h5ad

