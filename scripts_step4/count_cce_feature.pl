#!/usr/bin/perl -w
use strict;

my ($infile,$outfile) = @ARGV;

my %cce2feat;
my %cce2len;
open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$cce2feat{$array[3]}++;
	$cce2len{$array[3]}+=$array[8];
}
close IN;

open(OUT,">$outfile") || die $!;
print OUT "Clade\tFrequency\tLength\n";
foreach my $cce (sort keys %cce2feat) {
	print OUT $cce,"\t",$cce2feat{$cce},"\t",$cce2len{$cce},"\n";
}
close OUT;
