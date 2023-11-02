#!/user/bin/perl -w
use strict;
#####read gRNA#####
my $gRNA=$ARGV[3];
#my $gRNA="RetCand/data/mouse_mm10/refGene_RetGene_CandGene_crispr_target_1exonNotLast_100perc_screen_correct_kg_final.fa";
open(INPUT,$gRNA);
my %gRNA; ####gRNA hash
while(my $line = <INPUT>){ ####read one line
chomp $line;  
if($line =~ /^>(\S+)/){ ### if it is header line
my $line2 = <INPUT>; #####read the additional sequence line
chomp $line2;
$gRNA{"$1"}=$line2; #####the key is header line id
}
}

####read sequence .fa file
#my $NGS = "RetCand/data/*/Jiaxiong_Inj_1_S2_L001_R1_001.fa";
my $NGS = $ARGV[0]; ###############WT NGS fa file
my %NGS;
open(INPUT,$NGS);
while(my $line = <INPUT>){ ####read one line
chomp $line;
if($line =~ /^>(\S+)/){ ### if it is header line
my $line2 = <INPUT>; #####read the sequence line
chomp $line2;
$NGS{"$1"}=$line2; #### the key is header line id
}
}

#####read blastn result
#my $blastn = "RetCand/data/gRNA_pipeline/Jiaxiong_Inj_1_S2_L001_R1_001_gRNA";
my $blastn = $ARGV[1]; ######WT gRNA blastn result
open(INPUT,$blastn);
my %match;
my %NGS_match;
while(my $line = <INPUT>){ 
if($line =~ /^#/){ #####if start with #, skip the line
next;
}
my @info = split(/\s+/,$line); #### read each column
#if(($info[2] == 100) && ($info[3]==20)){ ##### if gRNA has 100% identity with sequencing and aligned length ==20, extract gRNA sequence and NGS sequence and check the match again
#if(($info[2] == 100) && ($info[3]>=14)){ ##### if gRNA has 100% identity with sequencing and aligned length ==14, extract gRNA sequence and NGS sequence and check the match again

#if($info[0] eq "emptyVector"){
#if(($info[2] == 100) && ($info[3]==37)){
#$match{$info[0]}->{$info[1]}=1;
#$NGS_match{$info[1]}=1;
#}
#}else{


#if(($info[2] == 100) && ($info[3]>=37)){ ##### if gRNA has 100% identity with sequencing and aligned length ==37, extract gRNA sequence and NGS sequence and check the match again

if(($info[2] == 100) && ($info[3]==26)){ 

#if(($info[2] >= 90) && ($info[3]>=65)){ ##### if gRNA has 100% identity with sequencing and aligned length ==37, extract gRNA sequence and NGS sequence and check the match again

#Abca4:chr3:+:122068965:122069162:6:50   @A00431:329:HGLL5DSX2:1:2678:26756:31923        100.000 20      0       0       1       20      79      98      0.002   40.1
#########my $gRNA_current = $gRNA{$info[0]}; ###gRNA ID is $info[0]
##########my $NGS_current = $NGS{$info[1]};  ####NGS ID in $info[1]
######if($NGS_current =~ /$gRNA_current/){
$match{$info[0]}->{$info[1]}=1;
$NGS_match{$info[1]}=1;
######}
}
#}
}

#my $output = $ARGV[2]; #####output WT gRNA count
my $output = $ARGV[2]."_new";
open(OUTPUT,">$output");
my $output1 = $ARGV[2]."_new_unmatch";
open(OUTPUT1,">$output1");
my $output2 = $ARGV[2]."_new_match";
open(OUTPUT2,">$output2");

#my $total_NGS = keys %NGS; #####total NGS reads
my $total_NGS_raw = keys %NGS;
my $total_NGS = keys %NGS_match;
for my $gRNA ( keys %gRNA){
if(defined $match{$gRNA}){
my $NGS_count = keys %{$match{$gRNA}};
my $NGS_frac = $NGS_count/$total_NGS;
my $NGS_frac_raw = $NGS_count/$total_NGS_raw;
print OUTPUT "$gRNA\t$NGS_count\t$NGS_frac\t$total_NGS\t$NGS_frac_raw\t$total_NGS_raw\n";
}else{
print OUTPUT "$gRNA\t0\t0\t$total_NGS\t0\t$total_NGS_raw\n";

#print OUTPUT "$gRNA\t0\t0\t$total_NGS\n";
}
}

for my $read ( keys %NGS){
if(!(defined $NGS_match{$read})){
print OUTPUT1 "$read\n$NGS{$read}\n";
}else{
print OUTPUT2 "$read\n";
}
}
