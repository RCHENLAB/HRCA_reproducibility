#!/usr/bin/perl -w
use strict;
my $dir="/storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL";
my $output1 = "$dir"."/all_TRUE_hg19_retina_gencodeHg19_1000g_flt_pip_total_merged_new";
my $output = "$dir"."/all_TRUE_hg19_retina_gencodeHg19_1000g_flt_cs_size_total_merged_new";
my $output2 = "$dir"."/all_TRUE_hg19_retina_gencodeHg19_1000g_flt_pip_var_merged_new";
open(OUTPUT,">$output");
open(OUTPUT1,">$output1");
open(OUTPUT2,">$output2");
my $output3 = "$dir"."/all_TRUE_hg19_retina_gencodeHg19_1000g_flt_cs_set_size_total_merged_new";
open(OUTPUT3,">$output3");
my $retnet = "cwc27/data/new_analysis/retnet_gene-1-7-2020_retcap_v5";
my %retnet;
open(INPUT,$retnet);
while(my $line = <INPUT>){
chomp $line;
$retnet{$line}=1;
}
my %cs_total=();
my %pip_total=();
my %cs_set_total=();
#my $anno="sc_human_retina/data/eQTL/all_TRUE_hg19_retina_gencodeHg19_anno_uniq";
#my $anno="sc_human_retina/data/cell_atlas/union_peak_bed_slope";
my $anno = "$dir"."/all_TRUE_hg19_retina_gencodeHg19.1000g_flt_snp_anno_new.gz";
open(INPUTa,"gunzip -c $anno | ");
my %anno;
<INPUTa>;
#10_100081757_C_A_b37    ENSG00000155229 10      100081757       8.05473e-05     -0.323568       0       9.3135974382319e-05     TRUE    10      99258551        ENSG00000155229 MMS19   NA      Non_peak        LD_block:chr10:97822357-100241302       -0.323568       99258551        823207  NotCRE  0.284869485949222       0.274644419630543       0.457712612939174
#       0.34606831822956        0.247150373718575       0.569623893502741       0.595709879259143       nonGene
while(my $linea = <INPUTa>){
chomp $linea;
my @infoa = split(/\s+/,$linea);
#RP11-455J20.3:2_190704420_G_A_b37       Non_peak        -0.583757
my @id = split(/\:/,$infoa[0]);
#my @id1 = split(/\_/,$id[1]);
my $snp = join(":", @id[0..5]); #"$id[0]:$id1[1]:$id1[2]:$id1[3]:$id[0]";
$anno{$snp} = $infoa[1];
}
#for(my $w=0; $w<=95;$w++){
#my $folder = "sc_human_retina/data/finemapping_eQTL/eQTL/all_TRUE_hg19_retina_gencodeHg19_gene100_1000g_flt_$w"."_other/";
my $folder = "sc_human_retina/data/finemapping_eQTL_new/all_TRUE_hg19_retina_gencodeHg19_other/";

####open(INPUT,"ls $folder/*prior_anno_sum | ");

open(INPUT,"ls $folder/*prior_uniform_sum | ");


while(my $line = <INPUT>){
my %cs=();

chomp $line;
my @info = split(/\//,$line);
my $cs = "$line";
my $gene = $info[-1];

$gene =~ s/\.prior_uniform_sum//;


####$gene =~ s/\.prior_anno_sum//;
#my $pip = "$folder/$info[-2].prior_pip";
my $pip = "$folder/$gene.prior_pip";

open(INPUT1,$pip) || die ("$pip");
my $header=<INPUT1>;
chomp $header;
my @header = split(/\s+/,$header);
my $id=0;
while(my $line1 = <INPUT1>){
chomp $line1;
#var     zscore  anno_pip        anno_lbf        uniform_pip     uniform_lbf
#chr9:76979797:rs79610609        -7.46   4.0956127378422e-13     25.2248604661636        3.33411076525181e-12    25.2248792148323
my @info1 = split(/\s+/,$line1);
$id++;
for(my $i=0; $i<=$#header; $i++){
$cs{$id}->{$header[$i]} = $info1[$i];
}
}
open(INPUT2,$cs);
<INPUT2>;
<INPUT2>;
<INPUT2>;
<INPUT2>;
print "$cs\n";
while(my $line2 = <INPUT2>){
chomp $line2;
print "$line2\n";
if((!($line2 =~ /\S+/)) || ($line2 =~ /^\<0/)){
last;
}
my @info2 = split(/\s+/,$line2);
my @var_id = split(/\:/,$cs{$info2[1]}->{var});
my $var_id = "$var_id[0]:$var_id[1]:$var_id[2]:$var_id[3]:$var_id[4]";
my $var_id_rev = "$var_id[0]:$var_id[1]:$var_id[3]:$var_id[2]:$var_id[4]"; 
$pip_total{$cs{$info2[1]}->{var}}= $cs{$info2[1]};
#print "$info2[1]\t$cs{$info2[1]}->{var}\n";
#exit;
#my $cs_name = "$info[-2]"."_$info2[-1]";
my $cs_name = "$gene"."_$info2[-1]";

$cs_total{$cs_name}->{$cs{$info2[1]}->{var}}=1;
#$cs_set_total{$info[-2]}->{$info2[-1]}=1;
$cs_set_total{$gene}->{$info2[-1]}=1;
print OUTPUT2 "$cs_name\t";
for(my $i=0; $i<=$#header; $i++){
print OUTPUT2 "$cs{$info2[1]}->{$header[$i]}\t";
}
if(defined $anno{$var_id}){
print OUTPUT2 "$anno{$var_id}\t";
}elsif(defined $anno{$var_id_rev}){
print OUTPUT2 "$anno{$var_id_rev}\t";
}
#if(defined $retnet{$var_id[4]}){
#print OUTPUT2 "YES";
#}else{
#print OUTPUT2 "NO";
#}
print OUTPUT2 "\n";
}
}
#}

my %pip_interval=();
my %cs_size=();
my %cs_set_size=();
for my $cs_set_key (keys %cs_set_total){
my $cs_set_num = keys %{$cs_set_total{$cs_set_key}};
$cs_set_size{$cs_set_num}++;
}

for my $cs_set_num ( keys %cs_set_size){
print OUTPUT3 "$cs_set_num\t$cs_set_size{$cs_set_num}\n";
}


for my $cs_key (keys %cs_total){
my $var_num = keys %{$cs_total{$cs_key}};
if($var_num==1){
$cs_size{size1}++;
}elsif(($var_num >=2) && ($var_num<=5)){
$cs_size{size2_5}++;
}elsif(($var_num >=6) && ($var_num<=10)){
$cs_size{size6_10}++;
}elsif(($var_num >=11) && ($var_num<=20)){
$cs_size{size11_20}++;
}elsif(($var_num >=21) && ($var_num<=50)){
$cs_size{size21_50}++;
}elsif(($var_num >=51)){
$cs_size{size51}++;
}
}

for my $var (keys %pip_total){
if(($pip_total{$var}->{anno_pip} >=0.01) && ($pip_total{$var}->{anno_pip} <0.1)){
$pip_interval{anno}->{perc1_10}++;
}

if(($pip_total{$var}->{uniform_pip} >=0.01) && ($pip_total{$var}->{uniform_pip} <0.1)){
$pip_interval{uniform}->{perc1_10}++;
}
###
if(($pip_total{$var}->{anno_pip} >=0.1) && ($pip_total{$var}->{anno_pip} <0.5)){
$pip_interval{anno}->{perc10_50}++;
}

if(($pip_total{$var}->{uniform_pip} >=0.1) && ($pip_total{$var}->{uniform_pip} <0.5)){
$pip_interval{uniform}->{perc10_50}++;
}
####
if(($pip_total{$var}->{anno_pip} >=0.5) && ($pip_total{$var}->{anno_pip} <0.8)){
$pip_interval{anno}->{perc50_80}++;
}

if(($pip_total{$var}->{uniform_pip} >=0.5) && ($pip_total{$var}->{uniform_pip} <0.8)){
$pip_interval{uniform}->{perc50_80}++;
}
#####
if(($pip_total{$var}->{anno_pip} >=0.8) && ($pip_total{$var}->{anno_pip} <0.9)){
$pip_interval{anno}->{perc80_90}++;
}

if(($pip_total{$var}->{uniform_pip} >=0.8) && ($pip_total{$var}->{uniform_pip} <0.9)){
$pip_interval{uniform}->{perc80_90}++;
}
#####
if(($pip_total{$var}->{anno_pip} >=0.9) && ($pip_total{$var}->{anno_pip} <0.95)){
$pip_interval{anno}->{perc90_95}++;
}

if(($pip_total{$var}->{uniform_pip} >=0.9) && ($pip_total{$var}->{uniform_pip} <0.95)){
$pip_interval{uniform}->{perc90_95}++;
}

####
if(($pip_total{$var}->{anno_pip} >=0.95) && ($pip_total{$var}->{anno_pip} <0.99)){
$pip_interval{anno}->{perc95_99}++;
}

if(($pip_total{$var}->{uniform_pip} >=0.95) && ($pip_total{$var}->{uniform_pip} <0.99)){
$pip_interval{uniform}->{perc95_99}++;
}
####
if(($pip_total{$var}->{anno_pip} >=0.99) ){
$pip_interval{anno}->{perc99}++;
}

if(($pip_total{$var}->{uniform_pip} >=0.99) ){
$pip_interval{uniform}->{perc99}++;
}
}

for my $key (keys %cs_size){
print OUTPUT "$key\t$cs_size{$key}\n";
}

for my $key (keys %{$pip_interval{anno}}){
print OUTPUT1 "anno\t$key\t$pip_interval{anno}->{$key}\n";
} 

for my $key (keys %{$pip_interval{uniform}}){
print OUTPUT1 "uniform\t$key\t$pip_interval{uniform}->{$key}\n";
}
