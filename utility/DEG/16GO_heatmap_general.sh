#!/bin/sh
err="human_meta/scripts/DEG/15GO_heatmap_general"
mkdir $err

for seq in all_new_all_mac_clean
do
for dem in gender ageGender age 
do
for ana in gsea #enricher
do
folder=${seq}_${dem}_${ana}
software/R-4.0.0/bin/Rscript --vanilla   DEG/16GO_heatmap_general.R $folder > ${err}/${folder}.out 2> ${err}/${folder}.err
done
done
done
