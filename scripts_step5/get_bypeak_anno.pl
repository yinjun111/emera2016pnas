#!/usr/bin/perl -w
use strict;

my %peaks;
open(IN,"cortex_ac_v5_all_bypeak_enh_phastConsElements46wayPlacental_all.txt") || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Enhancer/;
	my @array=split/\t/;
	$peaks{$array[0]}++;
}
close IN;

my %peak2anno;
my $title;
open(IN,"cortex_ac_v5_all_annotation_clade_all_forcount_140723.txt") || die $!;
while(<IN>) {
	tr/\r\n//d;
	if ($_=~/^Feature/) {
		$title=$_;
		next;
	}
	my @array=split/\t/;
	for(my $num=1;$num<@array;$num++) {
		$peak2anno{$array[0]}[$num-1]{$array[$num]}++;
	}
}
close IN;

open(OUT,">cortex_ac_v5_all_annotation_clade_all_forcount_140723_selbypeak.txt") || die $!;
print OUT $title,"\n";
foreach my $peak (sort keys %peaks) {
	my @contents;
	push @contents,$peak;
	foreach my $num (0..4) {
		push @contents,join(";",sort keys %{$peak2anno{$peak}[$num]});
	}
	print OUT join("\t",@contents),"\n";
}
close OUT;

		
