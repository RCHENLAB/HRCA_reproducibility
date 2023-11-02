#!/bin/sh
main_dir=/storage/chenlab/Users/junwang/monkey/scripts/database/
main_data_dir=/storage/chenlab/Users/junwang/monkey/data/database
export LD_LIBRARY_PATH=/storage/chen/Software/lib/usr/lib64:$LD_LIBRARY_PATH



file="/storage/chenlab/Users/junwang/human_meta/data/narrowPeaks_TimCherry_mac_ret_.bed"
file1="/storage/chenlab/Users/junwang/human_meta/data/narrowPeaks_TimCherry_mac_ret_hg19"

awk '{print $1":"$2"-"$3}' $file  > $file1

/storage/chen/Software/liftOver -positions  $file1  general_tool/hg19ToHg38.over.chain.gz  ${file1}_hg38 ${file1}_hg38_unmap


