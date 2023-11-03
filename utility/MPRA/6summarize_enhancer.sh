#!/bin/sh
#ll /storage/chenlab/Users/junwang/human_meta/data/proj2/h5ad/ | grep h5ad | awk '{print $NF}' | sed -e "s/_flt.h5ad//g" | sort -u > /storage/chenlab/Users/junwang/human_meta/data/proj2/snATAC_list
export PYTHONPATH=""
source /storage/chen/home/jw29/software/anaconda3/bin/activate  scvi-env
bc=$1
gRNA=$2
spe=$3
date=$4
python enhancer_validation/scripts/6summarize_enhancer.py $bc $gRNA $spe $date
