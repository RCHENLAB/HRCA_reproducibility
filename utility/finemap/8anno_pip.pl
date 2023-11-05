#!/usr/bin/perl 
use strict;
my $eQTL_file="sc_human_retina/data/eQTL/all_TRUE_hg19_retina_gencodeHg19";
my %eQTL=();
open(INPUT,$eQTL_file);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @coor=split(/\_/,$info[0]);
my $id="$coor[0]:$coor[1]:$coor[2]:$coor[3]"; 

$eQTL{$id}->{$info[-1]}=1;
$id="$coor[0]:$coor[1]:$coor[3]:$coor[2]"; 

$eQTL{$id}->{$info[-1]}=1;

}
my @cell = ("BC", "RGC", "Cone", "Rod", "MG", "HC","AC","Astrocyte","Microglia");
my %peak_all=();
my %peak_var=();
my $peak = $ARGV[0]; #####variant peak annotation
open(INPUT,$peak);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $peak = "chr$info[-4]:$info[-3]-$info[-2]";
$peak_all{$peak}=1;
my $pos = "$info[0]:$info[1]";
$peak_var{$pos}=$peak;
}

my %DAR;
for my $cell (@cell){
my $file = "/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/$cell"."_DAR_peak_hg38_hg19";
open(INPUT,$file);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
$DAR{$line}->{$cell}=1;
}
close(INPUT);
}


my %CRE;


my %CRE_gene;
#my $CRE_list = "/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/g2p_cA_final_major_full_max_LCRE_uniq_clean_hg38_hg19";
my $CRE_list = "/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/g2p_cA_final_major_full_max_LCRE_uniq_clean_hg38_hg19_out_CRE"; 
#hg38    hg19    gene
#chr10:100000226-100000726       chr10:101759983-101760483       DNMBP

open(INPUT,$CRE_list);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$CRE{$info[1]}->{$info[2]}=1;
}


my $peak_all = $ARGV[1]; ####variant annotation
my %peak;
my %gene;
open(INPUT,$peak_all);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\t/,$line);
$info[1] =~ s/chr//g;
my $var_pos = "$info[1]:$info[2]";
#print "$var_pos\n";
$gene{$var_pos} = $info[16];
if(($info[6] =~ /Intron/)||($info[6] =~ /Exon/)||($info[6] =~ /Promoter/)||($info[6] =~ /Downstream/)){
my @info1 = split(/\s+/,$info[6]);
$peak{$var_pos} = $info1[0];
}else{
my @info1 = split(/\s+/,$info[6]);
my $tmp = join("_",@info1);
$peak{$var_pos} = $tmp;
}
}

my $fn = $ARGV[2];
my $file = "/storage/chenlab/Users/junwang/human_meta/data/finemap/".$fn."_pip_var_all_merged_uniform";

#my $file = "/storage/chenlab/Users/junwang/human_meta/data/finemap/".$fn."_pip_var_all_merged";
my $output1=$file."_anno";
my $output2=$file."_finemap_gene";
my %finemap_gene;
my %finemap_gene_eQTL;
open(OUTPUT1,">$output1");
open(OUTPUT2,">$output2");
open(INPUT,$file);
while(my $line = <INPUT>){
#print "$line\n";
chomp $line;
#chr10:31907519-33707968_3       chr10:32094890:rs1023207:A:C    5.8033626454894 0.663494397659513       0.14458360653419   
my @info = split(/\s+/,$line);
my @coor = split(/\:/,$info[1]);

#######
#assign gene name
#######
my $var_pos = $coor[0].":".$coor[1];
#exit;
my $cate="";
my $assign=0;
$var_pos =~ s/chr//g ;
my $id=$coor[0].":".$coor[1].":".$coor[3].":".$coor[4];
#my $oppo_id = $coor[0].":".$coor[1].":".$coor[4].":".$coor[3];
#print "$var_pos\n";


if($peak{$var_pos} =~ /Exon/){
$cate.="Exon|"; ####label gene Exon
}

if($peak{$var_pos} =~ /3'_UTR/ ){
$cate.="3_UTR|"; #####label gene UTR
}

if($peak{$var_pos} =~ /5'_UTR/ ){
$cate.="5_UTR|"; #####label gene UTR
}


if($peak{$var_pos} =~ /Promoter/){
$cate.="Promoter|"; #####label gene UTR
}


if(exists $peak_var{$var_pos}){
my $peak = $peak_var{$var_pos};

if(exists $CRE{$peak}){
$cate.="CRE|"; #####label CRE
}


if(exists $DAR{$peak}){
$cate.="DAR|"; #####label CRE
}


$cate.="peak|"; #####label peak
}

my $gene = "";
my $eQTLL="";
$id =~ s/chr//g;
my %eQTL_gene;
my @gene;

if(exists $gene{$var_pos}){
$gene .= "var:$gene{$var_pos}|";
if($cate =~ /Exon|Promoter|3_UTR|5_UTR/){
push(@gene,$gene{$var_pos});

if(exists $eQTL{$id}){
if(exists $eQTL{$id}->{$gene{$var_pos}}){
$eQTLL.="$id:$gene{$var_pos}|";
$eQTL_gene{$gene{$var_pos}}=1;
}
}
}
}
if(exists $peak_var{$var_pos}){
my $peak = $peak_var{$var_pos};
if(exists $CRE{$peak}){
for my $key ( keys %{$CRE{$peak}}){
$gene .= "cre:$key|";
push(@gene,$key);
if(exists $eQTL{$id}){
if(exists $eQTL{$id}->{$key}){
$eQTLL.="$id:$key|";
$eQTL_gene{$key}=1;
}
}

}
}
}
print OUTPUT1 "$line\tanno;$cate\tgene;$gene\teQTL;$eQTLL\n";
#if((($cate =~ /[Exon|Promoter]/) && ($gene =~ /var/))|| (($cate =~ /CRE/) && ($gene =~ /cre/))){
#gene;var:LHX3|cre:QSOX2|cre:LHX3|
#my $gene1 = $gene;
#$gene1 =~ s/gene\;//g;
#$gene1 =~ s/var\://g;
#$gene1 =~ s/cre\://g;
#my @gene = split(/\|/,$gene1);
for my $gene (@gene){
$finemap_gene{$gene}=1;
if(exists $eQTL_gene{$gene}){
$finemap_gene_eQTL{$gene}=1;
}
}
#}

}
my $num=keys %finemap_gene;
my $num1=keys %finemap_gene_eQTL;
for my $gene ( keys %finemap_gene){
print OUTPUT2 "$gene\n";
}
print "gene:\t$num\ngene_eQTL:\t$num1\n";
 
