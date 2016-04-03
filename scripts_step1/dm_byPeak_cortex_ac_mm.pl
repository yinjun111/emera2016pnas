#!/usr/bin/perl -w
use strict;

#ind peaks
my %peaks;
open(IN,"../peaks/cortex_ac_mm_v5_all.bed") || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$peaks{$array[3]}++;
}
close IN;


#annotation of the bed files
open(IN,"bed_annotation_lcrv4sel_mm.txt") || die $!;
my %bed2anno1;
my %bed2anno2;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;	
	$bed2anno1{$array[0]}=$array[1];
	$bed2anno2{$array[0]}=$array[2];
}
close IN;

#original of the ind peaks from beds
my %peak2anno1;
my %peak2anno2;
foreach my $bedfile (sort keys %bed2anno1) {
	open(IN,"../peaks/$bedfile") || die $!;
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		my $bedindex=$array[0].":".$array[1]."-".$array[2];
		$peak2anno1{$bedindex}{$bed2anno1{$bedfile}}++;
		$peak2anno2{$bedindex}{$bed2anno2{$bedfile}}++;
	}
	close IN;
}
close IN;

open(OUT1,">cortex_ac_mm_v5_all_peakanno_140722.txt") || die $!;

my %merge2anno1;
my %merge2anno2;
print OUT1 "Feature\tByPeak1\tByPeak2\tByPeak3\n";
foreach my $peak (sort keys %peaks) {
	#anno by peak
	my %anno1;
	my %anno2;
	foreach my $a1 (sort keys %{$peak2anno1{$peak}}) {
		$anno1{$a1}++;
		#$merge2anno1{$ind2merge{$peak}}{$a1}++;
	}
	foreach my $a2 (sort keys %{$peak2anno2{$peak}}) {
		$anno2{$a2}++;
		#$merge2anno2{$ind2merge{$peak}}{$a2}++;
	}	
	
	#unless(defined $ind2merge{$peak}) {
	#	print STDERR $peak,"\n";
	#	next;
	#}
	
	#summary
	print OUT1 $peak,"\t";
	#ByPeak
	#By Peak
	if(keys %anno1==1) {
		print OUT1 (keys %anno1)[0],"\t";
	}
	else {
		print OUT1 "Multiple\t";
	}
	print OUT1 join(";",sort keys %anno1),"\t";	
	print OUT1 join(";",sort keys %anno2),"\n";	
}
close OUT1;
