#!/usr/bin/perl -w
use strict;

#my @files=glob("intersect/*_intersect.txt");

my ($infile,$outfile)=@ARGV;

my ($region, $clade);
if($infile=~/gencode.v10.wholeGene.exonTranscript_edi_([\dA-Za-z]+)_final_([A-Za-z]+)_intersect.txt/) {
	if($1 eq "10000") {
		$region="10k";
	}
	elsif($1 eq "100000") {
		$region="100k";
	}
	elsif($1 eq "100000s") {
		$region="100ks";
	}
	elsif($1 eq "intron") {
		$region="Intron";
	}	
	else {
		$region=$1;
	}
	$clade=$2;
	print STDERR $region,"\t",$clade,"\n";
}
else {
	print STDERR $infile,"\n";
	exit;
}

open(OUT,">$outfile") || die $!;
open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	#intersect cce(chr,s,e,name), region (chr,s,e,name)
	my ($istart,$iend,$ilen); #intersected cce
	if($array[1]<$array[5]) {
		$istart=$array[5];
	}
	else {
		$istart=$array[1];
	}
	if($array[2]>$array[6]) {
		$iend=$array[6];
	}
	else {
		$iend=$array[2]; #edited 14.1.6
	}
	$ilen=$iend-$istart;
	
	#the cce inside the region should be >=100bp, otherwise, it won't meet the conservation requirement
	if($ilen>=100) {
		print OUT join("\t",$array[0],$array[1],$array[2],$array[3],$region,$clade,$array[5],$array[6]),"\n";
	}
}
close IN;
close OUT;
#Homo_sapiens.GRCh37.67_genemerged_edi_mt100k_Amniota_intersect.txt