#!/usr/bin/perl -w
use strict;
my $gene_info = "sc_human_retina/data/eQTL/gencode.v19.genes.v7.patched_contigs.gtf_gene_tss_tes";
my $seq = $ARGV[0];
my %gene_info;
open(INPUTg,$gene_info);
while(my $lineg = <INPUTg>){
chomp $lineg;
my @infog = split(/\s+/,$lineg);
$gene_info{$infog[-1]}->{chr}=$infog[0];
$gene_info{$infog[-1]}->{start}=$infog[1];
$gene_info{$infog[-1]}->{end}=$infog[1]+1;
}
my %hash;

#my $sample = "/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_other_all2022_snRNA_rmHackney";
#my $sample = "/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_only_all_snRNA_young_mac_exp";
my $sample = "/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_other_all2023_mac_lobe_batch";
my %sample;
open(INPUTg,$sample);
my $lineg=<INPUTg>;
chomp $lineg;
my @sample_head=split(/\s+/,$lineg);
while(my $lineg = <INPUTg>){
chomp $lineg;
my @infog = split(/\s+/,$lineg);
for(my $i=0; $i<=$#infog; $i++){
$hash{$infog[2]}->{$sample_head[$i]}=$infog[$i];
}

#$hash{$infog[2]}->{$sample_head[$i]}=$infog[$i];
}
my @cell = ("AC","BC","Cone","HC","MG","RGC","Rod","Astrocyte"); #,"Microglia");

for my $cell (@cell){
#my $peer = "/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean/exp_".$cell."_logNorm_age_3peer";
#my $peer = "/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023/exp_".$cell."_logNorm_age_3peer";
my $peer = "/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_".$cell."_"."logNorm"."_".$seq."_3peer";

open(INPUT1,$peer);
<INPUT1>;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info = split(/\s+/,$line1);
my $label1 = "$cell"."_"."PC1";
my $label2 = "$cell"."_"."PC2";
my $label3 = "$cell"."_"."PC3";
$info[0] =~ s/Chen_b_D001.12/Chen_b_D001-12/g;
$info[0] =~ s/Roska_R.00646/Roska_R-00646/g;
$hash{$info[0]}->{$label1} = $info[1];
$hash{$info[0]}->{$label2} = $info[2];
$hash{$info[0]}->{$label3} = $info[3];


}
}
#Chr    start   end     TargetID        HG00096 HG00097 HG00099 HG00100 HG00101 HG00102 HG00103 HG00104 HG00105 HG00106 HG00108 HG00109 HG00110 HG00111 HG00112 HG00114 HG00115 HG00116 HG00117 HG00118 HG00119 HG00120
#D028_13 D027_13 D026_13 D021_13 D019_13 D018_13 D017_13 D013_13 D009_13 D005_13 D030_13 X19_D019        X19_D011        X19_D010        X19_D009        X19_D008        X19_D007        X19_D006        X19_D005        X19_D003        D19D015 D19D016 D19D013 D19D014
#A1BG    -1.7412906000252        0.180012369792705       1.27000841507014        0.627071687663682       -0.128239825980313      -0.449513565456691      0.901454079673677       -1.16283129020078       -0.0255806291232524     -0.565948821932863      2.04539098987329        0.393598415796924       -0.981125996162176      0.506871622557933       -1.39417320886912       1.54457548556325        0.757117297570861       -0.338888301553314      0.0768089757694082      -0.82713020584163       -0.690633454112706      -0.232272293482626      1.06757052387814        0.285174828345129

for my $cell (@cell){
#my $exp = "/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023/exp_".$cell."_logNorm_age";
my $exp = "/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_".$cell."_logNorm_".$seq;

my $output = "/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_".$cell."_logNorm_".$seq."_DEswan";

open(OUTPUT,">$output");
my @gene_array=();
open(INPUT,$exp);
my $header = <INPUT>;
$header =~ s/Chen_b_D001.12/Chen_b_D001-12/g;
$header =~ s/Roska_R.00646/Roska_R-00646/g;

chomp $header;
my @header = split(/\t/,$header);
my %gene_exp = ();

while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
for(my $i=1; $i<=$#info; $i++){

$gene_exp{"$header[$i-1]"}->{$info[0]} = $info[$i];
}
push(@gene_array, $info[0]);
}

my $label1 = "$cell"."_"."PC1";
my $label2 = "$cell"."_"."PC2";
my $label3 = "$cell"."_"."PC3";

my @var = ("$label1","batch","seq","age","race","gender");

#my @var = ("$label1","$label2","$label3","age","race","gender");
#print OUTPUT "PC1\tPC2\tPC3\tage\trace\tgender\t";
print OUTPUT "PC1\tbatch\tseq\tage\trace\tgender\t";

for(my $n=0; $n<=$#gene_array-1; $n++){
print OUTPUT "$gene_array[$n]\t";
}
print OUTPUT "$gene_array[$#gene_array]\n";

for my $ppl (sort keys %gene_exp){
print OUTPUT "$ppl\t";
$ppl =~ s/Chen_b_D001.12/Chen_b_D001-12/g;
$ppl =~ s/Roska_R.00646/Roska_R-00646/g;

for(my $j=0; $j<=($#var-1); $j++){
print OUTPUT "$hash{$ppl}->{$var[$j]}\t";

if(!(exists $hash{$ppl}->{$var[$j]}  )){
print "$cell\t$ppl\t$var[$j]\n";
}


}
if(!(exists $hash{$ppl}->{$var[$#var]}  )){
print "$cell\t$ppl\t$#var\n";
}
print OUTPUT "$hash{$ppl}->{$var[$#var]}";

for(my $n=0; $n<=$#gene_array; $n++){
print OUTPUT "\t$gene_exp{$ppl}->{$gene_array[$n]}";
if(!(defined $gene_exp{$ppl}->{$gene_array[$n]})){
print "$cell\tppl:$ppl\tgene:$gene_array[$n]\n";
exit;
}
}

print OUTPUT "\n";
}
}



