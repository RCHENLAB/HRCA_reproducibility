#!/usr/bin/perl -w
use strict;
my @cell = ("Rod","Cone","AC","BC","RGC","HC","MG","Astrocyte","Microglia");
my %cell_peak;
for my $cell (@cell){
#my $filename = "sc_human_retina/data/ATAC_peak/narrowPeak_combi/$cell"."_narrowPeak_combi_reform_macular_lobe_flt_1fpm";
my $filename="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/PeakCalls/".$cell."_peaks_census";
#my $filename = "/storage/chenlab/Users/junwang/sc_human_retina/data/snATAC_seq_peak/".$cell."_narrowPeak_combi_reform_macular_lobe_flt_2fpkm";
open(INPUT,$filename);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $peak = "$info[-4]:$info[-3]-$info[-2]";
if($info[-1] >0){
my $p = $info[-1]/($info[-2]-$info[-3]);
if($p >0.2){
$cell_peak{$peak}->{$cell}=1;
}
}
}
close(INPUT);
}

my %celltype_spe;
for my $peak ( keys %cell_peak){
my $cellnum = keys %{$cell_peak{$peak}};
$celltype_spe{$cellnum}->{"snATAC-seq OCRs"}->{$peak}=1;
}
my $output = "/storage/chenlab/Users/junwang/human_meta/data/peak_celltype";
open(OUTPUT,">$output");

my $overlap_bulk_ATAC= "/storage/chenlab/Users/junwang/human_meta/data/narrowPeaks_TimCherry_mac_ret_hg19_hg38_bed_cellatlas";

my %overlap_bulk_ATAC;
open(INPUT,$overlap_bulk_ATAC);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
#my $peak = "$info[0]:$info[1]-$info[2]";
if($info[-1] >0){
my $over_prop=$info[-1]/($info[-2]-$info[-3]);
if($over_prop > 0.2){
my $peak="$info[-4]:$info[-3]-$info[-2]";
#$overlap_bulk_ATAC{$peak}="$info[-4]:$info[-3]-$info[-2]";
my $cellnum = keys %{$cell_peak{$peak}};
$celltype_spe{$cellnum}->{"bulkATAC-seq OCRs"}->{$peak}=1;
}
}
}

open(OUTPUT,">$output");
for my $cellnum ( keys %celltype_spe){
for my $key ( keys %{$celltype_spe{$cellnum}}){
my $peak_num = keys %{$celltype_spe{$cellnum}->{$key}};
print OUTPUT "$cellnum\t$key\t$peak_num\n";
}
}
