#!/bin/sh
cell=$1
#cell="Both"
dir="human_meta/scripts/DEG/9pheatmap_age_snRNA_dream_batch_GO"
mkdir $dir


software/R-4.2.0/bin/Rscript --vanilla  human_meta/scripts/DEG/9pheatmap_age_snRNA_dream_batch_GO.R $cell >  ${dir}/${cell}.out 2> ${dir}/${cell}.err 
