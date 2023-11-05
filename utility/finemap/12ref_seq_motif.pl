#!/user/bin/perl -w
use strict;
#my $seq = "TAGACCGGACGCACGCCAAGCCGATGACCCACA";
#my $seq1 ="TAGACCGGACGCACGCTAAGCCGATGACCCACA";
#my $seq = uc("atccaTcaaattaattgccccaaatggga");
#my $seq1 = uc("atccaTcaaattaaCtgccccaaatggga");
#my $seq=uc("tcattttagaCaatgaaaggacattctttatttgctcct");
#my $seq = uc("ctCccagccatgacctcattgca");
my $seq = uc("gCcttcgtccggccgcagagg");
my $seq1 = uc("gCcttcgtccagccgcagagg");
#my $seq1 = uc("ctCccagccataacctcattgca");
#my $seq1=uc("tcattttagaCaatgaaagaacattctttatttgctcct");

#my $seq = uc("gctGacctggtcggtcatttaaatg");
#my $seq1 = uc("gctGacctggtcagtcatttaaatg");
#my $seq2=uc("tagaccggacacacgccaagccgatgacccaca");
#my $seq3=uc("tggaaccctgtgcccgccaagcctgatgaccaaca");
#my $seq4=uc("cggatcgcatgcagcccagcctgatgacccacg");
`mkdir sc_human_retina/scripts/ATAC_data/`;
my @nt = ("A","C","G","T");

#my $output = "sc_human_retina/scripts/ATAC_data/ref_PLEKHA7";
#my $output1 = "sc_human_retina/scripts/ATAC_data/alt_PLEKHA7";
#my $output = "sc_human_retina/scripts/ATAC_data/ref_HSF1";
#my $output1 = "sc_human_retina/scripts/ATAC_data/alt_HSF1";
#my $output = "sc_human_retina/scripts/ATAC_data/ref_RORA";
#my $output1 = "sc_human_retina/scripts/ATAC_data/alt_RORA";


my $output = "sc_human_retina/scripts/ATAC_data/ref_MBD2";
my $output1 = "sc_human_retina/scripts/ATAC_data/alt_MBD2";


#my $output = "sc_human_retina/scripts/ATAC_data/ref_NFE2L2";
#my $output1 = "sc_human_retina/scripts/ATAC_data/alt_NFE2L2";
#my $output2 = "sc_human_retina/scripts/ATAC_data/res_TTC9";
#my $output3 = "sc_human_retina/scripts/ATAC_data/mice_TTC9";
#my $output4 = "sc_human_retina/scripts/ATAC_data/dog_TTC9";
open(OUTPUT,">$output");
open(OUTPUT1,">$output1");

#open(OUTPUT2,">$output2");

#open(OUTPUT3,">$output3");

#open(OUTPUT4,">$output4");

#print OUTPUT "\t";
#print OUTPUT1 "\t";
for(my $i=1; $i<=length($seq); $i++){
print OUTPUT "\t$i";
print OUTPUT1 "\t$i";
#print OUTPUT2 "\t$i";
#print OUTPUT3 "\t$i";
#print OUTPUT4 "\t$i";

}
my @info = split(//,$seq);
my @info1 = split(//,$seq1);
print OUTPUT "\n";
print OUTPUT1 "\n";
for(my $i=0; $i<=$#nt; $i++){
print OUTPUT "$nt[$i]\t";
print OUTPUT1 "$nt[$i]\t";
for(my $j=0; $j<=$#info; $j++){
if($info[$j] eq $nt[$i]){
print OUTPUT "0.7\t";
}else{
print OUTPUT "0\t";
}

if($info1[$j] eq $nt[$i]){
print OUTPUT1 "0.7\t";
}else{
print OUTPUT1 "0\t";
}


}
print OUTPUT "\n";
print OUTPUT1 "\n";
}

exit;

#my @info2 = split(//,$seq2);
#for(my $i=1; $i<=length($seq2); $i++){
#print OUTPUT2 "\t$i";
#}
#print OUTPUT2 "\n";


#for(my $i=0; $i<=$#nt; $i++){
#print OUTPUT2 "$nt[$i]\t";
#for(my $j=1; $j<=$#info2; $j++){
#if($info2[$j] eq $nt[$i]){
#print OUTPUT2 "0.7\t";
#}else{
#print OUTPUT2 "0\t";
#}
#}
#print OUTPUT2 "\n";
#}
############3
#my @info3 = split(//,$seq3);
#for(my $i=1; $i<=length($seq3); $i++){
#print OUTPUT3 "\t$i";
#}
#print OUTPUT3 "\n";



#for(my $i=0; $i<=$#nt; $i++){
#print OUTPUT3 "$nt[$i]\t";
#for(my $j=1; $j<=$#info3; $j++){
#if($info3[$j] eq $nt[$i]){
#print OUTPUT3 "0.7\t";
#}else{
#print OUTPUT3 "0\t";
#}
#}
#print OUTPUT3 "\n";
#}

##########
#my @info4 = split(//,$seq4);
#for(my $i=1; $i<=length($seq4); $i++){
#print OUTPUT4 "\t$i";
#}
#print OUTPUT4 "\n";



#for(my $i=0; $i<=$#nt; $i++){
#print OUTPUT4 "$nt[$i]\t";
#for(my $j=1; $j<=$#info4; $j++){
#if($info4[$j] eq $nt[$i]){
#print OUTPUT4 "0.7\t";
#}else{
#print OUTPUT4 "0\t";
#}
#}
#print OUTPUT4 "\n";
#}



