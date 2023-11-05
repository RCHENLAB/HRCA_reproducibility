#!/usr/bin/perl -w
my $control_file = "/storage/chenlab/Users/junwang/human_meta/GWAS_file_list_new_control_MAGMA";
my %control;
open(INPUT,$control_file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$control{$info[0]} = $info[1];
}
my $output = $control_file."_summary_list_new";

open(OUTPUT,">$output");
print OUTPUT "P-val\tCell_type\tTrait\n";
my @list = ("MAGMA_list","MAGMA_list1");
for my $list (@list){
my $dir = "/storage/chen/home/jw29/human_meta/scripts/GWAS/$list";
open(INPUT,$dir);
while(my $line = <INPUT>){
chomp $line;

my $file = "$line/MAGMA_result_atlas"; #ad
open(INPUT1,$file)|| die ("$file");
<INPUT1>;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\t/,$line1);
$info1[1] =~ s/_MungeSumstats.txt.35UP.10DOWN//g;

$info1[1] =~ s/_MungeSumstats.txt_reform.35UP.10DOWN//g;
$info1[1] =~ s/_MungeSumstats.txt1_reform.35UP.10DOWN//g;
$info1[1] =~ s/_MungeSumstats_rmFRQ.txt_reform.35UP.10DOWN//g;
$info1[1] =~ s/_MungeSumstats.txt_rmFQS_reform.35UP.10DOWN//g;
if($info1[-3] eq "Linear"){
print OUTPUT "$info1[7]\t$info1[-1]\t$control{$info1[1]}\n";
}
}
}
}
 

