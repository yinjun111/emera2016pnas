#!/usr/bin/perl -w
use strict;

my %cate2color=("Human","230,230,250","Ape","230,230,250","Primate","230,230,250","Eutheria","230,230,250","Theria","230,230,250","Mammalia","255,0,0","Amniota","0,0,255","Tetrapoda","124,252,0","Gnathostomata","244,164,96","Vertebrata","0,0,0","NONE","192,192,192");

my @bedfiles=glob("merged/*.bed");


my %bed2info;
foreach my $bedfile (@bedfiles) {
	if($bedfile=~/hg19_nh.chrom.sizes_tobins100_50_v65_([A-Z\-\:a-z]+)_merged.bed/) {
		my $name=$1;
		my $cate;
		if($name=~/NONE/) {
			$cate="NONE";
		}
		else {
			$cate=$name;
		}
		open(IN,$bedfile) || die $!;
		while(<IN>) {
			tr/\r\n//d;
			my @array=split/\t/;
			my $chrindex=$array[0].":".$array[1]."-".$array[2];
			$bed2info{$chrindex}=[$array[0],$array[1],$array[2],$name,0,".",$array[1],$array[2],$cate2color{$cate}]
		}
		close IN;
	}
	else {
		print STDERR $bedfile,"\n";
	}
}

#color bed
open(OUT,">hg19_nh.chrom.sizes_tobins100_50_v65_merged_color.bed") || die $!;
foreach my $peak (sort {$bed2info{$a}[0] cmp $bed2info{$b}[0] || $bed2info{$a}[1] <=> $bed2info{$b}[1] } keys %bed2info) {
	print OUT join("\t",@{$bed2info{$peak}}),"\n";
}
close OUT;
