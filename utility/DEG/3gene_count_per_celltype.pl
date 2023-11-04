#!/usr/bin/perl -w
use strict;
my $rm_list = "/storage/chenlab/Users/junwang/human_meta/data/donor101_rm_new";
my %rm_list;
open(INPUT,$rm_list);
while(my $line = <INPUT>){
chomp $line;
$rm_list{$line}=1;
}
my $sample_list = "/storage/chenlab/Users/junwang/human_meta/data/donor_all_batch_new_snRNA_clean2023_all_celltype_num";


my %sample;
open(INPUT,$sample_list);
my $header = <INPUT>;
chomp $header;
my @sample_list1 = split(/\,/,$header);
my @sample_list = @sample_list1; #@sample_list1[0..30];
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\,/,$line);
for(my $i=1; $i<=$#info; $i++){
#if($info[$i]>=40){
if(defined $rm_list{$sample_list[$i]}){
next;
}
#if($info[$i]>=50){

if($info[$i]>=100){
$sample{$info[0]}->{$sample_list[$i]}=$info[$i];
}
}
}

my @ct_name = ("AC","BC","Cone","HC","MG","RGC","Rod","Astrocyte") ; #,"Microglia","RPE");
my %people;
for(my $s=1;$s<=$#sample_list;$s++){
my $sam = $sample_list[$s];
if(defined $rm_list{$sample_list[$s]}){
next;
}
#my $file = "/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_raw_batch_new_clean/$sam".".txt.gz";
my $file = "/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_raw_batch_new_clean_all/$sam".".txt.gz";

open(INPUT,"gunzip -c $file|");
my $line = <INPUT>;
chomp $line;
my @ct = split(/\s+/,$line);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
for(my $i=1; $i<=$#info;$i++){
$people{$ct[$i]}->{$info[0]}->{$sam} = $info[$i];
}
}
}
`mkdir /storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all`;

for my $key (keys %people){
if(($key eq "unassigned")) {#||($key eq "Rod")){
next;
}
my $output = "/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_"."$key";

open(OUTPUT,">$output");
my @sample_ct = sort (keys %{$sample{$key}});
my $line=join("\n", @sample_ct);
my $header = join("\t",@sample_ct);
print OUTPUT "$header\n";
for my $gene (sort keys %{$people{$key}}){
my $zero_gene=0;
    for(my $i=0; $i<=$#sample_ct; $i++){
        if(!(defined $people{$key}->{$gene}->{$sample_ct[$i]})){
                  $people{$key}->{$gene}->{$sample_ct[$i]}=0;
        $zero_gene++;                  
        }
    }

print OUTPUT "$gene\t";

    for(my $i=0; $i<=$#sample_ct-1; $i++){
        if(!(defined $people{$key}->{$gene}->{$sample_ct[$i]})){
                  $people{$key}->{$gene}->{$sample_ct[$i]}=0;
        }
     print OUTPUT "$people{$key}->{$gene}->{$sample_ct[$i]}\t";
    }


  my $tmp = $#sample_ct;
  if(!(defined $sample_ct[$#sample_ct])){
  print "$key\t$tmp\n";
  exit;
  }
  if(!(defined $people{$key}->{$gene}->{$sample_ct[$#sample_ct]})){
                  $people{$key}->{$gene}->{$sample_ct[$#sample_ct]}=0;
   }
     print OUTPUT "$people{$key}->{$gene}->{$sample_ct[$#sample_ct]}\n";
}
}
