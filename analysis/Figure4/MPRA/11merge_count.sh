#!/bin/sh
#ll /storage/chenlab/Users/junwang/human_meta/data/proj2/h5ad/ | grep h5ad | awk '{print $NF}' | sed -e "s/_flt.h5ad//g" | sort -u > /storage/chenlab/Users/junwang/human_meta/data/proj2/snATAC_list
export PYTHONPATH=""
source /storage/chen/home/jw29/software/anaconda3/bin/activate  scvi-env
python enhancer_validation/scripts/11merge_count.py
