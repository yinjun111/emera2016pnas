#!/usr/bin/perl -w
use strict;
use List::Util qw(min max);

my ($infile,$outfile)=@ARGV;
#infile, enhancer intersect with cce, wo
#outfile, bed file of the intersection

my %names;
open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	unless(defined $names{$array[3].";".$array[7]}) {
		print OUT $array[0],"\t",max($array[1],$array[5]),"\t",min($array[2],$array[6]),"\t",$array[3],";",$array[7],"\n";
		$names{$array[3].";".$array[7]}++;
	}
	else {
		print OUT $array[0],"\t",max($array[1],$array[5]),"\t",min($array[2],$array[6]),"\t",$array[3],";",$array[7],";",$names{$array[3].";".$array[7]},"\n";
		$names{$array[3].";".$array[7]}++;
	}		
}
close IN;
close OUT;
