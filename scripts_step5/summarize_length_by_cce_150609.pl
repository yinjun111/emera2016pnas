#!/usr/bin/perl -w
use strict;
use List::Util qw(min max);

#output a matrix for each enhancer
#clade;eutheria->mtamn with best LOD score


#summarzie the conservation score for enhancers

#cortex_ac_v5_all_bypeak_enh_140723.bed
#phylop.vertebrate


#read in enh annotation #do calculation for CS16, CS23, F2 separately
my %clade2contents;
my %clade2name=qw(Human Human Ape Ape Primate Primate Amniota Amniota Gnathostomata MtAmn Mammalia Mammalia Tetrapoda MtAmn Vertebrata MtAmn Eutheria Eutheria Theria Theria MtAmn MtAmn NONE:Amniota NONE NONE:Tetrapoda NONE NONE:Gnathostomata NONE NONE:Vertebrata NONE NONE NONE EutheriaNOlder EutheriaNOlder TheriaNOlder TheriaNOlder MammaliaNOlder MammaliaNOlder AmiotaNOlder AmniotaNOlder);
my @names=qw(Human Ape Primate Eutheria Theria Mammalia Amniota MtAmn);
my @namessel=qw(EutheriaNOlder TheriaNOlder MammaliaNOlder AmniotaNOlder MtAmn);


my %peaks;
open(IN,"cortex_ac_v5_all_annotation_clade_all_forcount_140723_selbypeak.txt") || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Feature/;
	my @array=split/\t/;
	$peaks{$array[0]}++;
}
close IN;

#cce to len
my %enh2cce;
my $enh2cce_file="cortex_ac_v5_all_bypeak_enh_140723_intersect_cceaggremerged.bed";
open(IN,$enh2cce_file) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	my ($enh,$clade)=split(";",$array[3]);
	if(defined $clade2name{$clade}) {
		$enh2cce{$enh}{$clade2name{$clade}}+= $array[2]-$array[1];
	}
	else {
		print STDERR $clade,"\n";
	}
}
close IN;


#enh to len
my %enh2len;
open(IN,"cortex_ac_v5_all_bypeak_enh_140723.bed") || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$enh2len{$array[3]}= $array[2]-$array[1];
}
close IN;


print STDERR "Here!\n";


my $outfile="cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_length_bycceaggresummary.txt";
open(OUT,">$outfile") || die $!;
print OUT "Enhancer\tTotalLength\t",join("\t", @namessel),"\n";
foreach my $enh (sort keys %peaks) {
	print OUT $enh,"\t";
	if(defined $enh2len{$enh}) {
		print OUT $enh2len{$enh},"\t";
	}
	else {
		if($enh=~/chr\w+\:(\d+)-(\d+)/) {
			print OUT $2-$1,"\t";
		}
	}
	my @scores;
	foreach my $clade (@namessel) {
		if(defined $enh2cce{$enh}{$clade}) {
			push @scores,$enh2cce{$enh}{$clade};
		}
		else {
			push @scores,0;
		}
	}
	print OUT join("\t",@scores),"\n";
}
close OUT;
