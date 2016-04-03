#!/usr/bin/perl -w
use strict;

my ($infile,$bedfile,$outfile)=@ARGV;

my %clade2peak;
open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Peak/;
	my @array=split/\t/;
	#unless($array[4] =~/NONE/) {
		$clade2peak{$array[4]}{$array[0]}++;
	#}
	#else {
	#	$clade2peak{"NONE"}{$array[0]}++;
	#}
}
close IN;

my %bed2coord;
open(IN,$bedfile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	push @{$bed2coord{$array[3]}},join("\t",@array[0..2]);	
}
close IN;

foreach my $clade (sort keys %clade2peak) {
	open(OUT,">$outfile\_$clade.bed") || die $!;
	foreach my $bed (sort keys %{$clade2peak{$clade}}) {
		foreach my $coord (@{$bed2coord{$bed}}) {
			print OUT join("\t",$coord,$bed,$clade),"\n";  #there will be some regions overlapping with each other with different conservation level
		}
	}
	close OUT;
}
