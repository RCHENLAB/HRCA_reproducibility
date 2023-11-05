#!/usr/bin/perl -w

my @peak = ("peak");
for my $peak (@peak){
my $dir = "/storage/chenlab/Users/junwang/human_meta/data/GWAS/major"."_$peak/ldsc_result/";

#my $dir = "/storage/chenlab/Users/junwang/human_meta/data/GWAS/group"."_$peak/ldsc_result/";
my $list = "human_meta/scripts/GWAS/file_list_new_control"; ####final peak


my $output = $list."_$peak"."_summary";
my $output1 = $list."_$peak"."_summary_list";

open(OUTPUT,">$output");
open(OUTPUT1,">$output1");
print OUTPUT1 "P-val\tCell_type\tTrait\tCoef\n";
my %label;
my %hash;
open(INPUT,$list);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
#my @info1 = split(/\//,$info[0]);
my $file = "$dir/$info[0]".".cell_type_results.txt";
if(open(INPUT1,$file)){ # || next ; # die "$file";
$label{$info[1]}=1;

<INPUT1>;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info2 = split(/\s+/,$line1);
$hash{$info2[0]}->{$info[1]}=$info2[3];
print OUTPUT1 "$info2[3]\t$info2[0]\t$info[1]\t",$info2[1]*1000000,"\n";
}
}
}
my @key = keys %label;
for my $key (sort @key){
print OUTPUT "$key\t";
}
print OUTPUT "\n";
for my $key (sort keys %hash){
print OUTPUT "$key\t";
for my $key1 (sort @key){
print OUTPUT -log($hash{$key}->{$key1})/log(10),"\t";
#print OUTPUT "$hash{$key}->{$key1}\t";
}
print OUTPUT "\n";
}
}
