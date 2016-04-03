#!/usr/bin/perl -w
use strict;

my ($infile,$bedfile,$outfile)=@ARGV;

my @testlins=qw(Human Ape Primate Eutheria Theria Mammalia Amniota Tetrapoda Gnathostomata Vertebrata);
my @testspecs=qw(hg19	panTro2	ponAbe2	rheMac2	calJac1	mm9	rn4	cavPor3	oryCun2	bosTau4	equCab2	canFam2	loxAfr3	monDom5	ornAna1	galGal3	taeGut1	anoCar1	xenTro2	tetNig2	fr2	gasAcu1	oryLat2	danRer6	petMar1);
my @matchedlins=qw(Human Ape Ape Primate Primate Eutheria Eutheria Eutheria Eutheria Eutheria Eutheria Eutheria Eutheria Theria Mammalia Amniota Amniota Amniota Tetrapoda Gnathostomata Gnathostomata Gnathostomata Gnathostomata Gnathostomata Vertebrata);

#my %clade2col=(Human=>[1],Ape=>[2,3],Primate=>[4,5],Eutheria=>[6..13],Theria=>[14],Mammalia=>[15],Amniota=>[16..18],Tetrapoda=>[19],Gnathostomata=>[20..24],Vertebrata=>[25]);
    
#record all the peaks, so peaks without alignment will be assigned as human
#peak_forSubmit
my %peaks;
open(IN,$bedfile) || die $!;
while(<IN>) {
	tr/\r\n//d;
    my @array=split/\t/;
	$peaks{$array[3]}++;
}
close IN;

#read processed age results
open(IN,$infile) || die $!;
my %peak2info;

while(<IN>) {
    tr/\r\n//d;
    my @array=split/\t/;
    next if $_=~/^Peak/;
    $peak2info{$array[0]}=$_;
}
close IN;

#print oUT results
open(OUT,">$outfile") || die $!;
print OUT "Peak\tCategory\t","MostDistantlyConsistentlyID\tMostDistantlyConsistentlySpec\tMostDistantlyConsistentlyClade\t","MostDistantlyID\tMostDistantlySpec\tMostDistantlyClade\t";
print OUT join("\t",@testspecs),"\t",join("\t",@testlins),"\n";
foreach my $peak (sort keys %peaks) {
	if(defined $peak2info{$peak}) {
		print OUT $peak2info{$peak},"\n";
	}
	else {
		print OUT $peak,"\tNONE\t1\thg19\tHuman\t1\thg19\tHuman\t";
		print OUT join("\t",1,(0) x (scalar(@testspecs)-1)),"\t",join("\t",1,(0) x (scalar(@testlins)-1)),"\n";
	}
}
close OUT;



