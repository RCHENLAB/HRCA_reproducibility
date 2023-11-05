#!/usr/bin/perl
use strict;
#my $dir = "sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_sig_snp_1000g_flt_20ppl";
my $dir = $ARGV[0];
my $dir1 = $dir."_other";
`rm -r $dir1`;
`mkdir $dir1`;
open(INPUT, "ls $dir | ");
#my $zscore="sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8.1000g_flt_20ppl.zscore.gz";
my $zscore = "$ARGV[0].zscore.gz";
my %zscore;
open(INPUT1,"gunzip -c $zscore |");
while(my $line1 = <INPUT1>){
chomp $line1;
#chr6:29027255:rs112383057       chr6:28917608-29737971  -5.977
my @info1 = split(/\s+/,$line1);
$zscore{$info1[0]} = $info1[2];
}
while(my $line = <INPUT>){
chomp $line;
#print "line:$line\n";
#exit;
my $file="$dir/$line";
my $sort ="$dir1/$line"."_sort";
my $snp = "$dir1/$line"."_snp";
my @info = split(/\//,$line);
my @coor = split(/\:/,$info[-1]);
my $chr  = $coor[0];
$chr =~ s/chr//g;
#generate a tmp file
`less $file | sed -e "s/\\:/\\t/g" | sort -k 1,1 -k 2n,2n > $sort`;
`awk '{print \$3}' $sort > $snp`;
#print "$sort\n$snp\n";
#exit;
my $plink="/storage/chen/Software/plink_linux_x86_64/plink";
my $input="/storage/chen/home/jw29/sc_human_retina/data/GWAS/1000G_EUR_Phase3_plink/1000G.EUR.QC."."$chr";
my $out="$dir1/$line"."_r";
`$plink --bfile $input --extract $snp --r square spaces  --keep-allele-order  --out $out`;
my $zout = "$dir1/$line"."_zscore";
open(OUTPUTz,">$zout");
open(INPUTs,"$sort");
while(my $lines = <INPUTs>){
chomp $lines;
my @infos = split(/\s+/,$lines);
my $id = join(":",@infos[0..4]);
my $zscore = $zscore{$id};
print OUTPUTz "$id\t$zscore\n";
}
}
#`cut -d ":" $file | sort -k 1,1 -k 2,2 > $tmp`;
