#!/usr/bin/perl -w
use strict;

my ($infile,$annofile,$outfile,$region)=@ARGV;
#region can be "intron", or a number

my %gene2exon;
open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$gene2exon{$array[3]}{$array[1]}=[@array[0..2]];
}
close IN;

my %chr2size;
open(IN,$annofile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$chr2size{$array[0]}=$array[1];
}
close IN;


open(OUT,">$outfile") || die $!; #intron
foreach my $gene (sort keys %gene2exon) {
	my @starts=sort {$a <=> $b } keys %{$gene2exon{$gene}};
	#intron
	if($region eq "intron") {
		if(@starts>1) {
			for(my $num=1;$num<@starts;$num++) {
				if($gene2exon{$gene}{$starts[$num-1]}[2]<$gene2exon{$gene}{$starts[$num]}[1]) {
					print OUT $gene2exon{$gene}{$starts[0]}[0],"\t",$gene2exon{$gene}{$starts[$num-1]}[2],"\t",$gene2exon{$gene}{$starts[$num]}[1],"\t",$gene,"\n";
				}
			}
		}
	}
	elsif($region=~/^\d+$/) {
		#regions near the gene
		my $start=$gene2exon{$gene}{$starts[0]}[1];
		my $end=$gene2exon{$gene}{$starts[$#starts]}[2];
		my ($bstart,$bend,$cstart,$cend);
		
		if($start>=$region) {
			$bstart=$start-$region;
			$bend=$start;
		}
		else {
			$bstart=0;
			$bend=$start;
		}
		if($bstart<$bend && $bend>0) {
			print OUT $gene2exon{$gene}{$starts[0]}[0],"\t",$bstart,"\t",$bend,"\t",$gene,"\n";
		}

		if($end+$region<=$chr2size{$gene2exon{$gene}{$starts[0]}[0]}) {
			$cstart=$end;
			$cend=$end+$region;
		}
		else {
			$cstart=$end;
			$cend=$chr2size{$gene2exon{$gene}{$starts[0]}[0]};
		}
		if($cstart<$cend && $cend>0) {
			print OUT $gene2exon{$gene}{$starts[0]}[0],"\t",$cstart,"\t",$cend,"\t",$gene,"\n";
		}
	}
}
close OUT;
