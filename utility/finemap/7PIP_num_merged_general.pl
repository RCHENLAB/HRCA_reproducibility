#!/usr/bin/perl -w
use strict;
my $folder = "$ARGV[0]";
my $dir="/storage/chenlab/Users/junwang/human_meta/data/finemap";
`mkdir $dir`;
my $output1 = "$dir/$ARGV[1]"."_pip_total_merged";
my $output = "$dir/$ARGV[1]"."_cs_size_merged";
my $output2 = "$dir/$ARGV[1]"."_pip_var_all_merged";


open(OUTPUT,">$output");
open(OUTPUT1,">$output1");
open(OUTPUT2,">$output2");
open(INPUT,"ls $folder/*prior_anno_sum | ");
my %cs_total;
my %pip_total;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\/|\./,$line);
my $cs = "$line";
#print "$cs";
#exit;
my $pip = "$folder/$info[-2].prior_pip";
my %cs=();
open(INPUT1,$pip);
my $header=<INPUT1>;
chomp $header;
my @header = split(/\s+/,$header);
my $id=0;
while(my $line1 = <INPUT1>){
chomp $line1;
#var     zscore  anno_pip        anno_lbf        uniform_pip     uniform_lbf
#chr9:76979797:rs79610609        -7.46   4.0956127378422e-13     25.2248604661636        3.33411076525181e-12    25.2248792148323
#var     zscore  anno_pip        anno_lbf        uniform_pip     uniform_lbf
#chr9:76207149:rs9314808:A:T     5.83509477532426        0.000171087733782604    14.7307630258651        0.000206803570145797    14.7293115752898
#chr9:76207300:rs9314809:A:G     5.84351726019672        0.000179532212427169    14.778941186125 0.000217012934555116    14.7774991270522
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
while(my $line2 = <INPUT2>){
chomp $line2;
print "$line2\n";
if((!($line2 =~ /\S+/)) || ($line2 =~ /^\<0/)){
last;
}
my @info2 = split(/\s+/,$line2);
# variable variable_prob cs
#       30     0.9999329  1
#print "0:$info2[0]\t1:$info2[1]\n";
#exit;
$pip_total{$cs{$info2[1]}->{var}}= $cs{$info2[1]};
my $cs_name = "$info[-2]"."_$info2[-1]";
$cs_total{$cs_name}->{$cs{$info2[1]}->{var}}=1;
print OUTPUT2 "$cs_name\t";

for(my $i=0; $i<=$#header; $i++){
print OUTPUT2 "$cs{$info2[1]}->{$header[$i]}\t";
}
print OUTPUT2 "\n";
}
}

my %pip_interval;
my %cs_size;

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
