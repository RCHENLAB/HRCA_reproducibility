#!/usr/bin/perl -w
use strict;
my @cell =("Rod", "BC", "AC", "Cone", "HC", "MG", "RGC");
my %hash;
my %hash_up;
my %hash_down;
my $dem="age";
my $seq="cpm07_all_new_all_mac_clean";
for my $cell (@cell){

my $file = "/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/$cell"."_interval_cpm01_snRNA_clean_young_".$seq."/".$cell."_DEG_res_cpm_".$dem;

open(INPUT,$file);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
if(($info[-1] ne "NA" ) &&($info[-1] < 0.05)){

$hash{$info[0]}->{$cell}=1;
if($info[1] > 0){

$hash_up{$info[0]}->{$cell}=1;
}

if($info[1] < 0){

$hash_down{$info[0]}->{$cell}=1;
}

}



}
}
my %sum;
for my $gene ( keys %hash){
my @cell1= sort keys %{$hash{$gene}};
my $cell_line = join("&",@cell1);
$sum{$cell_line}->{$gene}=1;
}

my %sum_up;
for my $gene ( keys %hash_up){
my @cell1= sort keys %{$hash_up{$gene}};
my $cell_line = join("&",@cell1);
$sum_up{$cell_line}->{$gene}=1;
}


my %sum_down;
for my $gene ( keys %hash_down){
my @cell1= sort keys %{$hash_down{$gene}};
my $cell_line = join("&",@cell1);
$sum_down{$cell_line}->{$gene}=1;
}

my $dir1 = "/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream";
my $output = "$dir1/upset_list_$dem"."_all_".$seq;

my %count;

open(OUTPUT,">$output");
my $line="";
for my $key ( keys %sum){
my $gene_num = keys %{$sum{$key}};
$line .= "\'$key\'=$gene_num,";
my @key=split("&",$key);
if($#key==0){
$count{$key}->{only}=$gene_num;
$count{$key}->{all}+=$gene_num;
}else{
for my $key1 ( @key){
$count{$key1}->{all}+=$gene_num;
}
}
}


my $new_line = substr($line,0,length($line)-1);
print OUTPUT "$new_line";

my $output10="$dir1/upset_list_$dem"."_all_".$seq."_ct";
open(OUTPUT,">$output10");
for my $key (keys %count){
my $prop = $count{$key}->{only}/$count{$key}->{all};
print OUTPUT "$key\t$count{$key}->{only}\t$count{$key}->{all}\t$prop\n";
}

my $output = "$dir1/upset_list_$dem"."_up_".$seq;

open(OUTPUT,">$output");
my $line="";
for my $key ( keys %sum_up){
my $gene_num = keys %{$sum_up{$key}};
$line .= "\'$key\'=$gene_num,";
}

my $new_line = substr($line,0,length($line)-1);
print OUTPUT "$new_line";

my $output = "$dir1/upset_list_$dem"."_down_".$seq;
open(OUTPUT,">$output");
my $line="";
for my $key ( keys %sum_down){
my $gene_num = keys %{$sum_down{$key}};

$line .= "\'$key\'=$gene_num,";
}

my $new_line = substr($line,0,length($line)-1);
print OUTPUT "$new_line";
