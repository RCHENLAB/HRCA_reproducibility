#!/usr/bin/perl
use strict;
my $main = "/storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL";
my $dir = "$main"."/all_TRUE_hg19_retina_gencodeHg19";
my $dir1 = $dir."_other";
`rm -r $dir1`;
`mkdir $dir1`;
open(INPUT, "ls $dir | ");
my $zscore = "$main"."/all_TRUE_hg19_retina_gencodeHg19.1000g_flt_matrixeQTL_formt_new.gz"; 
my %zscore;
open(INPUT1,"gunzip -c $zscore |");
<INPUT1>;
while(my $line1 = <INPUT1>){
chomp $line1;
#chr6:29027255:rs112383057       chr6:28917608-29737971  -5.977
my @info1 = split(/\s+/,$line1);
$zscore{$info1[0]} = $info1[3];
}
while(my $line = <INPUT>){
chomp $line;
my $file="$dir/$line";
my $sort ="$dir1/$line"."_sort";
my $snp = "$dir1/$line"."_snp";
my $tmp = "$dir1/$line"."_tmp";
my @info = split(/\//,$line);

`less $file | sed -e "s/\\:/\\t/g" | sort -u  | sort -k 1,1 -k 2n,2n > $sort`;

`cut -f 1 $sort | sort -u > $tmp`;
`awk '{print \$6}' $sort > $snp`;

open(INPUTt,$tmp);
my $chr  = <INPUTt>;
chomp $chr;
my $plink="/storage/chen/Software/plink_linux_x86_64/plink";
my $input="/storage/chen/home/jw29/sc_human_retina/data/GWAS/1000G_EUR_Phase3_plink/1000G.EUR.QC."."$chr";
my $out="$dir1/$line"."_r";
`$plink --bfile $input --extract $snp --r square  --out $out`;
my $zout = "$dir1/$line"."_zscore";
open(OUTPUTz,">$zout");
open(INPUTs,"$sort");
while(my $lines = <INPUTs>){
chomp $lines;
my @infos = split(/\s+/,$lines);
#my $id = join(":",@infos[0..2]);
my $id = join(":",@infos[0..5]);
my $zscore = $zscore{$id};
print OUTPUTz "$id\t$zscore\n";
}
}
#`cut -d ":" $file | sort -k 1,1 -k 2,2 > $tmp`;
