#!/usr/bin/perl -w
use strict;
my $file = "sc_human_retina/data/eQTL/all_TRUE_hg19_retina_gencodeHg19";
my $vcf = $file.".vcf";
#my $tmp_out = "$vcf"."_out";
#my $tmp="sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/all_peak_anno_bin_bed";
my $tmp ="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/peakset_clusters_major_full_peak_anno_hg38_hg19_bed";
`echo "##fileformat=VCFv4.2" > $vcf`;
`echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG001">> $vcf`;
#`awk '{split(\$1,b,\"\_\")\; split(\$2,a,\"\:\")\; print \"chr\"a[1]\"\\t\"a[2]\"\\t\"a[5]\"\\t\"a[3]\"\\t\"a[4]\"\\t\.\\tPASS\\t\"b[1]\"\\tGT\\t1/1\"}' sc_human_retina/data/finemapping_eQTL/eQTL/all_TRUE_hg19_retina_gencodeHg19_gene100_1000g_flt_pip_var_merged >> $vcf`;
`tail -n+2 $file | awk '{split(\$1,a,\"\_\")\; print a[1]\"\\t\"a[2]\"\\t\.\t\"a[3]\"\\t\"a[4]\"\\t\.\\tPASS\\t\"\$NF\"\\tGT\\t1/1\"}' | sort -k 1,1 -k 2n,2n  >> $vcf`;


my $tmp_out = $file."_format_eQTL_atlas";
`bedtools intersect -a $vcf -b $tmp -wo > $tmp_out`;

