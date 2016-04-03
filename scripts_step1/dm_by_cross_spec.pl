#!/usr/bin/perl -w
use strict;

my ($hsdiff,$mmdiff,$reci,$result,$comparecol,$spec) = @ARGV;

unless(defined $comparecol) {
	$comparecol=2; #real column number
}
unless(defined $spec) {
	$spec="inter"; #or unique #unique tissue, or multiple tissue allowed
}


my %human2anno;
my @title1;
open(IN,$hsdiff) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	if($_=~/^Feature/) {
		@title1=@array[1..$#array];
	}
	else {
		$human2anno{$array[0]}=[@array[1..$#array]];
	}
}
close IN;

my %mouse2anno;
my @title2;
open(IN,$mmdiff) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	if($_=~/^Feature/) {
		@title2=@array[1..$#array];
	}
	else {
		$mouse2anno{$array[0]}=[@array[1..$#array]];
	}
}
close IN;

my %human2mouse;
open(IN,$reci) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	next unless defined $mouse2anno{$array[1]}; #sometime the reci annotation is redundant
	$human2mouse{$array[0]}{$array[1]}++;
}
close IN;

open(OUT,">$result") || die $!;
print OUT "Feature\tByBoth\tByHuman\tByHuman\tByMouse\tByMouse\tMousePeaks\t",join("\t",@title1),"\t",join("\t",@title2),"\n";
foreach my $peak (sort keys %human2anno) {
	print OUT $peak,"\t";
	if(defined $human2mouse{$peak}) {
		my %bymouse;
		my @mouseannos;
		foreach my $mpeak (sort keys %{$human2mouse{$peak}}) {
			unless($spec eq "unique") {
				foreach my $anno (split(";",$mouse2anno{$mpeak}[$comparecol-2])) {
					$bymouse{$anno}++;
				}
			}
			else {
				$bymouse{$mouse2anno{$mpeak}[$comparecol-2]}++;
			}
			for(my $num=0;$num<scalar(@{$mouse2anno{$mpeak}});$num++) {
				push @{$mouseannos[$num]},$mouse2anno{$mpeak}[$num];
			}
		}
		my %byhuman;
		unless($spec eq "unique") {
			foreach my $anno (split(";",$human2anno{$peak}[$comparecol-2])) {
				$byhuman{$anno}++;
			}
		}
		else {
			$byhuman{$human2anno{$peak}[$comparecol-2]}++;
		}		
		
		#if human anno exist in mouse anno, then it is byBoth
		my %byboth;		
		foreach my $anno (sort keys %byhuman) {
			if(defined $bymouse{$anno} ) {
				$byboth{$anno}++;
			}
			elsif($anno eq "hESC" && defined $bymouse{"mESC"}) {
				$byboth{$anno}++;
			}
			elsif($anno eq "small_intestine" && defined $bymouse{"intestine"}) {
				$byboth{$anno}++;
			}
		}
		
		if(keys %byboth>0) {
			print OUT join(";",sort keys %byboth),"\t";
		}
		else {
			print OUT "Nomatch\t";
		}
			
		#human anno
		if(keys %byhuman>1) {
			print OUT "Multiple\t";
		}
		else {
			print OUT (keys %byhuman)[0],"\t";
		}
		
		print OUT join(";",sort keys %byhuman),"\t";
		
		#mouse anno
		if(keys %bymouse==1) {
			print OUT (keys %bymouse)[0],"\t";
		}
		else {
			print OUT "Multiple\t";
		}
		print OUT join(";",sort keys %bymouse),"\t";
		
		#mouse peaks
		print OUT join(";",sort keys %{$human2mouse{$peak}}),"\t";
		#human all
		print OUT join("\t",@{$human2anno{$peak}}),"\t";
		#mouse all
		print OUT join("\t",map {join(";",@$_)} @mouseannos),"\n";
	}
	else {
		#ByBoth
		print OUT "None\t";
		#ByHuman
		#human anno
		my %byhuman;
		unless($spec eq "unique") {
			foreach my $anno (split(";",$human2anno{$peak}[$comparecol-2])) {
				$byhuman{$anno}++;
			}
		}
		else {
			$byhuman{$human2anno{$peak}[$comparecol-2]}++;
		}		
		if(keys %byhuman>1) {
			print OUT "Multiple\t";
		}
		else {
			print OUT (keys %byhuman)[0],"\t";
		}
		
		print OUT join(";",sort keys %byhuman),"\t";
		
		#mouse anno
		print OUT "None\t";
		print OUT "None\t";
				
		#mouse peaks
		print OUT "None\t";
		#human all
		print OUT join("\t",@{$human2anno{$peak}}),"\t";
		#mouse all
		print OUT join("\t", ("None") x @{$human2anno{$peak}}),"\n";
	}
}
close OUT;




