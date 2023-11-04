#!/bin/sh
export PYTHONPATH=""

for seq in all_new_all_mac_clean
do
for dem in age gender ageGender
do

folder="${seq}_${dem}"
dem="${folder}_gsea"
ct=1000000
DEG/15GO_gsea_general.py  $folder $dem $ct

done
done

