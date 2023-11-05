#!/usr/bin/perl 
use strict;
use Math::CDF qw(qnorm);
my %G1000_var;

my %G1000_var;
for(my $i=1; $i<=22; $i++){
my $bim="sc_human_retina/data/eQTL/all_TRUE_hg19_retina_gencodeHg19_chr".$i."_1000G.bim";
#9       rs4741076       25.284641       10978031        T       G
#1       rs55727773      0.37021441      706368  A       G
open(INPUTb,$bim);
while(my $lineb = <INPUTb>){
chomp $lineb;
my @info = split(/\s+/,$lineb);
my $var = "$info[0]:$info[3]:$info[4]:$info[5]";
$G1000_var{$var}=$info[1];
}
}

my @cell = ("BC", "RGC", "Cone", "Rod", "MG", "HC","AC","Astrocyte","Microglia");
my %peak_all=();
my %peak_var=();
my $peak = "sc_human_retina/data/eQTL/all_TRUE_hg19_retina_gencodeHg19_format_eQTL_atlas"; #####variant peak annotation
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
#my $file = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/$cell"."_DAR_07-23-2021_flt_bin_new";
my $file = "/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/$cell"."_DAR_peak_hg38_hg19";
open(INPUT,$file);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
#my @info = split(/\t/,$line);
#my $peak = "$info[0]:$info[1]-$info[2]";
#$DAR{$peak}->{$cell}=1;
$DAR{$line}=1;
}
close(INPUT);
}
#my $peak_all_anno ="/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/all_peak_anno_bin_new";

#my %peak_gene;
#my %peak_anno;
#open(INPUT,$peak_all_anno);
#<INPUT>;
#while(my $line = <INPUT>){
#chomp $line;
#my @info = split(/\t/,$line);

#my $var_pos = "$info[1]:$info[2]-$info[3]";
#$info[1] =~ s/chr//g;
#$peak_gene{$var_pos} = $info[16];
#if(($info[6] =~ /Intron/)||($info[6] =~ /Exon/)||($info[6] =~ /Promoter/)||($info[6] =~ /Downstream/)){
#my @info1 = split(/\s+/,$info[6]);
#$peak_anno{$var_pos} = $info1[0];
#}else{
#my @info1 = split(/\s+/,$info[6]);
#my $tmp = join("_",@info1);
#$peak_anno{$var_pos} = $tmp;
#}
#}


my %CRE;

my %CRE_gene;
#my $CRE_list = "sc_human_retina/data/snATAC_snRNA/LCRE_list";
my $CRE_list = "/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/g2p_cA_final_major_full_max_LCRE_uniq_clean_hg38_hg19";
open(INPUT,$CRE_list);
while(my $line = <INPUT>){
chomp $line;
#my @info = split(/\s+/,$line);
#my $peak = "chr$info[1]:$info[2]-$info[3]";
#my $link = $peak."_".$info[0];
$CRE{$line}=1;
#$CRE_gene{$info[0]}->{$peak}=1;
}


my $peak_all = "sc_human_retina/data/eQTL/all_TRUE_hg19_retina_gencodeHg19_format_anno"; ####variant annotation
my %peak;
my %gene;
open(INPUT,$peak_all);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\t/,$line);
$info[1] =~ s/chr//g;
my $var_pos = "$info[1]:$info[2]";
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


my $file = "sc_human_retina/data/eQTL/all_TRUE_hg19_retina_gencodeHg19";
`mkdir /storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL/`;
my $out="/storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL/all_TRUE_hg19_retina_gencodeHg19";
my $output1=$out.".1000g_flt_matrixeQTL_formt_new";
my $output2=$out.".1000g_flt_snp_anno_new";
my $output3=$out.".1000g_flt_snp_map_new";
my $output4=$out.".1000g_flt_gene_map_new";
open(OUTPUT1,">$output1");
open(OUTPUT2,">$output2");
open(OUTPUT3,">$output3");
open(OUTPUT4,">$output4");
print OUTPUT1 "SNP\tgene\tbeta\tt-stat\tp-value\n";
print OUTPUT2 "SNP\tannot_d\n";



#chr:6:29027255  a       g       0.4933  0.0068  232207.00       5.977   2.271e-09       ++++    0.0     1.056   3       0.7876
open(INPUT,$file);
<INPUT>;
#<INPUT>;
while(my $line = <INPUT>){
chomp $line;
#variant_id	gene_id	variant_chr	variant_pos	pval	slope	flag	genelevel_threshold	is_signif	gene_chr	gene_tss	gene_id_add	genename
#10_1584564_T_C_b37	ENSG00000151240	10	1584564	3.30854e-05	0.336801	1	7.75984046504918e-05	TRUE	10	735683	ENSG00000151240	DIP2C
my @info = split(/\s+/,$line);
my $sign_dir;
if($info[5] > 0){
$sign_dir =1;
}elsif($info[5]<0){
$sign_dir=-1;
}


my $egene = $info[-1];
my $egene_tss = $info[-3];
my @coor = split(/_/,$info[0]);
my $zscore = $sign_dir*abs(qnorm($info[4]/2));
my $beta = $info[5];
my $chr = "chr".$coor[0]; #"$coor[0]"."$coor[1]";
my $var = "$coor[0]:$coor[1]:".uc($coor[2]).":".uc($coor[3]);
my $oppo_var = "$coor[0]:$coor[1]:".uc($coor[3]).":".uc($coor[2]);



my $rs;
my $snp;
if(defined $G1000_var{$var}){
$rs=$G1000_var{$var};
$snp=$var;
}elsif(defined $G1000_var{$oppo_var}){
$zscore = -$zscore;
$beta = -$beta;
$rs=$G1000_var{$oppo_var};
$snp=$oppo_var;
}else{
next;
}


print OUTPUT1 "$snp:$info[-1]:$rs\t$info[-1]\t$beta\t$zscore\t$info[4]\n";

print OUTPUT3 "$snp:$info[-1]:$rs\t$coor[0]\t$coor[1]\n";

print OUTPUT4 "$egene\t$coor[0]\t$egene_tss\n"; 

#`rm $tmp`;
#######
#assign gene name
#######
my $var_pos = $coor[0].":".$coor[1];
my $cate=0;
my $assign=0;

if(($peak{$var_pos} =~ /Exon/)&&($gene{$var_pos} eq $info[-1])){
$cate=5; ####label gene Exon
$assign=1;
}

if(($assign==0)&&(($peak{$var_pos} =~ /3'_UTR/ )||( $peak{$var_pos} =~ /5'_UTR/  ))&&( $gene{$var_pos} eq $info[-1]  )){
$cate=5; #####label gene UTR
$assign=1;
}


if(($assign==0)&&($peak{$var_pos} =~ /Promoter/)&&( $gene{$var_pos} eq $info[-1]  )){
$cate=4; #####label gene UTR
$assign=1;
}

if(defined $peak_var{$var_pos}){
my $peak = $peak_var{$var_pos};


if(($assign==0)&& (defined $CRE{$peak})){ # && (defined $CRE{$peak}->{$info[-1]})){
$cate=3; #####label CRE
$assign=1;
}

if(($assign==0)&& (defined $DAR{$peak})){ # && (defined $CRE{$peak}->{$info[-1]})){
$cate=2; #####label CRE
$assign=1;
}

if($assign==0){
$cate=1; #####label peak
$assign=1;
}
}

if($assign==0){
$cate=0;
$assign=1;
}

print OUTPUT2 "$snp:$info[-1]:$rs\t$cate\n";


}

`bgzip $output1`;
`bgzip $output2`;
`bgzip $output3`;
`bgzip $output4`;

