#!/bin/sh
for seq in all_new_all_mac_clean
do
for cell in Cone BC Astrocyte HC MG RGC Rod AC
do
dem="ageGender"

sbatch -p interactive --mem=5000MB /storage/chen/home/jw29/human_meta/scripts/DEG/9pheatmap_age_snRNA_dream_batch_GO_general.sh $cell $seq $dem

dem="age"

sbatch -p interactive --mem=5000MB /storage/chen/home/jw29/human_meta/scripts/DEG/9pheatmap_age_snRNA_dream_batch_GO_general.sh $cell $seq $dem

dem="gender"

sbatch -p interactive --mem=5000MB /storage/chen/home/jw29/human_meta/scripts/DEG/9pheatmap_age_snRNA_dream_batch_GO_general.sh $cell $seq $dem


done
done
