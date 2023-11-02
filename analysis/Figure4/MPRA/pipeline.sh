#!/bin/sh

date="5-17-23"
file="file_list_hm3"
spe="hm"
dir_tmp="stat1"

sh 0file_list.sh $date

sh 1format_query.sh $date 

sh 2format_gRNA_database.sh $date

sh 4blastn_enhancer_final_hm.sh  $date

sh 5count_perfect_match_gRNA_hm.sh $date

sh 6read_count_dist_12422_log.sh $date $dir_tmp

####merge count#######
sh 11merge_count.sh

####summarize dropoff#######
sh 6summarize_enhancer_batch.sh  $date $file $spe

#####process count######
sh 8process_count.sh $date $file

######barcode mean#######
sh 9barcode_mean.sh $date $spe

######analyze#######
sh 10compair_basal.sh $date $spe $dir_tmp

######plot 
sh 12plot_IRD.R
