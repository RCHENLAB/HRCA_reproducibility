#!/usr/bin/perl -w 
use strict;
my $main = "/storage/chen/home/jw29";
#my @file =("sc_human_retina/data/GWAS/SummaryStat/coloc_retina_eQTL_sum_stringent-9-14-2023_PPH4_0.5_pip_VSX2_RHO"); 
#my @file = ("sc_human_retina/data/GWAS/SummaryStat/coloc_retina_eQTL_sum_stringent-9-14-2023_PPH4_0.5_pip_GLCCI1");
#my @file = ("sc_human_retina/data/GWAS/SummaryStat/coloc_retina_eQTL_sum_stringent-9-14-2023_PPH4_0.5_pip_GLCCI1");
my @file = ("sc_human_retina/data/GWAS/SummaryStat/coloc_retina_eQTL_sum_stringent-10-27-2021_PPH4_0.5_HTRA1_CRE_hg19_coloc");
#my @file = ("sc_human_retina/data/GWAS/SummaryStat/coloc_retina_eQTL_sum_stringent-10-27-2021_PPH4_0.5_pip_poag");
#my @file = ("sc_human_retina/data/GWAS/SummaryStat/coloc_retina_eQTL_sum_stringent-10-27-2021_PPH4_0.5_overlap_QTL_ASCA");
#my $dir_zoomlocus= "sc_human_retina/data/GWAS_coloca/Zoomlocus_retina";
#`mkdir $dir_zoomlocus`;
my $gene_file = "sc_human_retina/data/eQTL/gencode.v19.genes.v7.patched_contigs.gtf_gene_tss_tes";
my %gene_name;
my %list;
open(INPUT,$gene_file);
#1       11869   14362   ENSG00000223972.4       DDX11L1
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @ensembl = split(/\./,$info[-2]);
$gene_name{$info[-1]}=$ensembl[0];
}
for my $file(@file){
open(INPUT,$file);
#sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB/          sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB//coloc_retina_bulk_eQTL//10_111835691_ENSG00000119950_summary    
while(my $line = <INPUT>){
chomp $line;
my $dir_zoomlocus;
my $pre_loc_file;
my @info = split(/\t/,$line);
my $dir= $info[0];
#my @var = split(/\//,$line); #####to get the variant name 
#my @info1 = split(/\_/,$var[-1]);
my @info1 = split(/\_/,$info[1]);
#my $rs = $info[24];
my $chr = $info1[0];
my $pos = $info1[1];
my $rs = "$chr"."_".$pos;
my $gene = $info1[2];
if($dir =~ /Khawaja2018_29785010/){
$pre_loc_file = "sc_human_retina/data/GWAS/SummaryStat/Khawaja2018_29785010/Khawaja2018_29785010NG/UKBBc.IOP.txt.rmDup.sumstat.P.sort.gz";
$dir_zoomlocus= "sc_human_retina/data/GWAS/SummaryStat/Khawaja2018_29785010/Khawaja2018_29785010NG/locusZoom_retina_coloc";

}elsif( $dir =~ /FritscheLG2016_26691988/){
$pre_loc_file = "sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform.sumstat.P.sort.gz";
$dir_zoomlocus= "sc_human_retina/data/GWAS/locusZoom_retina_coloc";
}elsif($dir =~ /Myopia29808027/){
$pre_loc_file = "/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/Myopia29808027/Myopia29808027NG/41588_2018_127_MOESM14_ESM.sorted.P.gz";
$dir_zoomlocus= "sc_human_retina/data/GWAS/SummaryStat/Myopia29808027/Myopia29808027NG/locusZoom_retina_coloc";
}elsif($dir =~ /KhorCC2016_27064256/){
$pre_loc_file = "sc_human_retina/data/GWAS/SummaryStat/KhorCC2016_27064256/KhorCC2016_27064256NG/Khor-27064256_reform.sumstat.P.sort.gz";
$dir_zoomlocus= "sc_human_retina/data/GWAS/SummaryStat/KhorCC2016_27064256/KhorCC2016_27064256NG/locusZoom_retina_coloc";
}elsif($dir =~ /ChoquetH_29891935/){
#print "enter\n";
$pre_loc_file = "sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB/POAG_GERA_UKB_meta_reform.sumstat.P.sort.gz";
$dir_zoomlocus= "sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB/locusZoom_retina_coloc";
}elsif($dir =~ /CraigJ_31959993/){
$pre_loc_file = "sc_human_retina/data/GWAS/SummaryStat/CraigJ_31959993/CraigJ_31959993NG/MTAG_glaucoma_four_traits_summary_statistics.txt.sumstat.P.sort.gz";
$dir_zoomlocus= "sc_human_retina/data/GWAS/SummaryStat/CraigJ_31959993/CraigJ_31959993NG/locusZoom_retina_coloc";
}elsif($dir =~ /Hysi_32231278/){
$pre_loc_file = "sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020.txt.gz_rs_all_pos_sort_P.gz";
$dir_zoomlocus="sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/locusZoom_retina_coloc";
}elsif($dir =~ /Gharahkhani_33627673/){
$pre_loc_file = "sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37.sumstat.P.sort.gz";
$dir_zoomlocus = "sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/locusZoom_retina_coloc";
}elsif($dir =~ /36848389PG_ONL/){
$pre_loc_file = "sc_human_retina/data/GWAS/SummaryStat/ONL/Currant_36848389PG_ONL/GCST90243953_BuildGRCh37.tsv.P.sort.gz";
$dir_zoomlocus = "sc_human_retina/data/GWAS/SummaryStat/ONL/Currant_36848389PG_ONL/locusZoom_retina_coloc";
}elsif($dir =~ /36848389PG_OST/){
$pre_loc_file = "sc_human_retina/data/GWAS/SummaryStat/OST/Currant_36848389PG_OST/GCST90243955_BuildGRCh37.tsv.P.sort.gz";
$dir_zoomlocus = "sc_human_retina/data/GWAS/SummaryStat/OST/Currant_36848389PG_OST/locusZoom_retina_coloc";
}


if($pre_loc_file eq "NA"){
next;
}
#my $pos_start = $pos -500000;
#my $pos_end = $pos + 500000;
#my $pos_start = $pos -200000;
#my $pos_end = $pos + 200000;
my $pos_start = $pos -50000;
my $pos_end = $pos + 50000;

my $temp_rs = "$main/$dir_zoomlocus/tmp_rs";
my $region = "chr$chr:$pos-$pos";

`tabix $pre_loc_file $region | cut -f 3 > $temp_rs`;
if(open(INPUTrs,$temp_rs)){
$rs = <INPUTrs>;
chomp $rs;
}
print "$rs\n";
`mkdir -p $dir_zoomlocus/$rs`;

my $gwas = "$main/$dir_zoomlocus/$rs/$rs"."_GWAS";
print "$gwas\n";
my $prefix = "$rs"."_GWAS_EUR500kb";
$list{"$gwas\t$chr\t$pos_start\t$pos_end\t$prefix\t$rs\t$dir_zoomlocus"}=1;
my $region = "chr$chr:$pos_start-$pos_end";
`echo "CHR	POS	SNP	A1	A2	P" > $gwas`;
`tabix $pre_loc_file $region | sed -e 's/chr//g' | cut -f 1-6 >> $gwas `;


my $retina_eQTL = "sc_human_retina/data/eQTL/Retina.nominal.eQTLs.with_thresholds/all_hg19_rs.sort.gz";
#V1      V2      V3      V4      V5      V6      V7      V8      V9      V10     V11     Effect_allele   Baseline_allele V12     V13     V14     genelevel_threshold     is_signif
#ENSG00000227232 chr1    29571   29570   -       983     19394   1:10177:A:AC    chr1    10177   10177   AC      A       0.381708
#        -0.0839295      0       0.000386075983656031    FALSE
my $eQTL ="$main/$dir_zoomlocus/$rs/$rs"."_$gene"."_eQTL";
print "$eQTL\n";

#my $prefix = "$dir_zoomlocus/$rs/$rs"."_$genename"."_eQTL_EUR200kb";
$prefix = "$rs"."_$gene"."_eQTL_EUR500kb";
$list{"$eQTL\t$chr\t$pos_start\t$pos_end\t$prefix\t$rs\t$dir_zoomlocus"}=1;
`echo "CHR      POS     SNP     A1      A2      P" > $eQTL`;
$region =~ s/chr//g;
`tabix $retina_eQTL $region | grep $gene | cut -f 1-6 >> $eQTL`;
}
}

#my $outputlist = "/storage/chen/home/jw29/sc_human_retina/data/GWAS_coloca/locusZoom_list_coloc_retina-5-13-2021";
my $outputlist = "/storage/chen/home/jw29/sc_human_retina/data/GWAS_coloca/locusZoom_list_coloc_retina-9-22-2022";

#my $outputlist = "/storage/chen/home/jw29/sc_human_retina/data/GWAS_coloca/locusZoom_list_coloc_retina-4-10-2021";
open(OUTPUTl,">$outputlist");
for my $key (keys %list){
print OUTPUTl "$key\n";
my @key = split(/\s+/,$key);
# `/storage/chen/Software/miniconda2/bin/python   /storage/chen/Software/locuszoom/bin/locuszoom --plotonly --verbose --metal sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB_6_36592987-36592987  --chr  6  --start 36342987 --end 36842987 --prefix r61105472EUR250kb  --build hg19 --pop EUR --source 1000G_Nov2014 --delim whitespace --pvalcol P --markercol SNP --snpset NULL --epacts-chr-col CHR --refsnp rs61105472`;
my $tmp_dir = "$key[6]/$key[5]";
 `cd $tmp_dir; /storage/chen/Software/miniconda2/bin/python   /storage/chen/Software/locuszoom/bin/locuszoom --plotonly --verbose --metal $key[0]  --chr  $key[1]  --start $key[2] --end $key[3] --prefix $key[4]  --build hg19 --pop EUR --source 1000G_Nov2014 --delim whitespace --pvalcol P --markercol SNP --snpset NULL  --refsnp $key[5] --epacts-chr-col CHR `;

}
