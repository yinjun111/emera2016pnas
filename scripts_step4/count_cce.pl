#!/usr/bin/perl -w
use strict;

my $infile="/Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/MamCortex/genome/Human/cce_v65/hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmnMerged.bed";
my $outfile="ccemerged_summary.txt";

my %cce2feat;
my %cce2len;
open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$cce2feat{$array[3]}++;
	$cce2len{$array[3]}+=$array[2]-$array[1];
}
close IN;

open(OUT,">$outfile") || die $!;
print OUT "Clade\tFrequency\tLength\n";
foreach my $cce (sort keys %cce2feat) {
	print OUT $cce,"\t",$cce2feat{$cce},"\t",$cce2len{$cce},"\n";
}
close OUT;
