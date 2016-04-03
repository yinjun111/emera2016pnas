#!/usr/bin/perl -w
use strict;

my %cce2anno;
my %cce2bklen;
my %lin2background;
open(IN,"gencode_v10_cce_all_shuffle_background_short.txt") || die $!;
while(<IN>){
	tr/\r\n//d;
	my @array=split/\t/;
	$cce2anno{$array[5]}{$array[0]}{$array[4]}{$array[3]}=[@array];
	$lin2background{$array[5]}{$array[0]}{$array[4]}{$array[3]}=$array[7]-$array[6]; #lin, chr,reg, name => len
	$cce2bklen{$array[3]}=$array[7]-$array[6];
}
close IN;


#assign index to regions
open(OUT,">gencode_v10_cce_all_shuffle_background_forShuffle_140606.txt") || die $!;
my %background2index;
foreach my $lin (sort keys %lin2background) {
	foreach my $chr (sort keys %{$lin2background{$lin}}) {
		foreach my $reg (sort keys %{$lin2background{$lin}{$chr}}) {
			my $indexnum=1;
			foreach my $cce (sort {$lin2background{$lin}{$chr}{$reg}{$b} <=> $lin2background{$lin}{$chr}{$reg}{$a}} keys %{$lin2background{$lin}{$chr}{$reg}}) {
				#$background2index{$lin}{$chr}{$reg}{$indexnum}=$cce;
				print OUT join("\t",@{$cce2anno{$lin}{$chr}{$reg}{$cce}},$indexnum),"\n";
				$indexnum++;
			}
		}
	}
}
close OUT;