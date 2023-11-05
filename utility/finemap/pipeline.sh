#!/bin/sh
file=$1
file_name=$2
sample_size=$3
#1.finemap
dir="utility/finemap"
export PERL5LIB=/storage/chen/home/jw29/software/perl/lib/perl5/Math-CDF-0.1/lib/perl5/x86_64-linux/:$PERL5LIB

#err=${dir}/1prepare_torus_profile

sh	${dir}/1prepare_torus_profile_perl.sh  ${file} ${file}_format_atlas_peak ${file}_anno

#2 run torus
command=${dir}/2run_torus_all.sh
#dir="human_meta/scripts/finemap/2run_torus_all"
#mkdir $dir

sh	$command ${file}_format_atlas.vcf.1fpm.1000g_flt_20ppl


#3
command=${dir}/3generate_r_plink_general

sh	${command}.sh ${file}_format_atlas.vcf.1fpm.1000g_flt_20ppl


#4
command=${dir}/4suise_general
sh	${command}.sh ${file}_format_atlas.vcf.1fpm.1000g_flt_20ppl $sample_size


#5
command=${dir}/5pip_plot_general_text

sh	${command}.sh ${file}_format_atlas.vcf.1fpm.1000g_flt_20ppl

#6

err=${dir}/7PIP_num_merged_general
command=${err}_perl

sh ${command}.sh ${file}_format_atlas.vcf.1fpm.1000g_flt_20ppl_other  ${file_name}

