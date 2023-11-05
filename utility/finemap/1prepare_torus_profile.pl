#!/usr/bin/perl 
use strict;
use Math::CDF qw(qnorm);
my %G1000_var;
for(my $i=1; $i<=22; $i++){
my $bim = $ARGV[0]."_chr$i"."_1000G.bim";
#9       rs4741076       25.284641       10978031        T       G
open(INPUTb,$bim);
while(my $lineb = <INPUTb>){
chomp $lineb;
my @info = split(/\s+/,$lineb);
my $var = "chr$info[0]:$info[3]:$info[4]:$info[5]";
$G1000_var{$var}=$info[1];
}
}

my @cell = ("BC", "RGC", "Cone", "Rod", "MG", "HC","AC","Astrocyte","Microglia");
my %peak_all=();
my %peak_var=();
my $peak = $ARGV[1]; #####variant peak annotation
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
my $CRE_list = "/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/g2p_cA_final_major_full_max_LCRE_uniq_clean_hg38_hg19";
open(INPUT,$CRE_list);
while(my $line = <INPUT>){
chomp $line;
$CRE{$line}=1;
}


my $peak_all = $ARGV[2]; ####variant annotation
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


my $file = "$ARGV[0]_format_atlas.vcf";
my $ldetect = "sc_human_retina/data/GWAS_new/ldetect/EUR/fourier_ls-all.bed.gz"; 
my $output1=$file.".1fpm.1000g_flt_20ppl.zscore";
my $output2=$file.".1fpm.1000g_flt_20ppl.anno";

my $output1_gz = "$output1".".gz";
my $output2_gz = "$output2".".gz";
`rm $output1_gz`;
`rm $output2_gz`;
open(OUTPUT1,">$output1");
open(OUTPUT2,">$output2");
print OUTPUT2 "SNP\tannot_d\n";
#chr:6:29027255  a       g       0.4933  0.0068  232207.00       5.977   2.271e-09       ++++    0.0     1.056   3       0.7876
open(INPUT,$file);
<INPUT>;
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $sign_dir;
if($info[5] > 0){
$sign_dir =1;
}elsif($info[5]<0){
$sign_dir=-1;
}

my $zscore = $sign_dir*abs(qnorm($info[6]/2));

my $chr = "chr$info[0]";
my $var = "$chr:$info[1]:".uc($info[4]).":".uc($info[3]);
my $oppo_var = "$chr:$info[1]:".uc($info[3]).":".uc($info[4]);

my $rs;
if(defined $G1000_var{$var}){
$rs=$G1000_var{$var};
}elsif(defined $G1000_var{$oppo_var}){
$zscore = -$zscore;
$rs=$G1000_var{$oppo_var};
}else{
next;
}
my $tmp = $ARGV[0]."_tmp";
`rm $tmp`;
`less sc_human_retina/data/GWAS_new/ldetect/EUR/fourier_ls-all.bed.gz | awk '{if((\$1==\"$chr\")&&(\$2<=$info[1])&&(\$3>=$info[1])){print \$1\":\"\$2\"-\"\$3}}' > $tmp`;
open(INPUT1,"$tmp");
my $ld = <INPUT1>;
chomp $ld;
print OUTPUT1 "$chr:$info[1]:$rs:$info[4]:$info[3]\t$ld\t$zscore\n";
#`rm $tmp`;
#######
#assign gene name
#######
my $var_pos = $info[0].":".$info[1];
my $cate=0;
my $assign=0;

if($peak{$var_pos} =~ /Exon/){
$cate=5; ####label gene Exon
$assign=1;
}

if(($assign==0)&&(($peak{$var_pos} =~ /3'_UTR/ )||( $peak{$var_pos} =~ /5'_UTR/  ))){
$cate=5; #####label gene UTR
$assign=1;
}
if(($assign==0)&&($peak{$var_pos} =~ /Promoter/)){
$cate=4; #####label gene UTR
$assign=1;
}


if(defined $peak_var{$var_pos}){
my $peak = $peak_var{$var_pos};

if(($assign==0)&& (defined $CRE{$peak})){
$cate=3; #####label CRE
$assign=1;
}


if(($assign==0)&& (defined $DAR{$peak})){
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

print OUTPUT2 "$chr:$info[1]:$rs:$info[4]:$info[3]\t$cate\n";




}

`bgzip $output1`;
`bgzip $output2`;

