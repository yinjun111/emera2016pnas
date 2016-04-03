#!/usr/bin/perl -w
use strict;
use List::Util qw(min max);

#summarzie the conservation score for enhancers

my ($bedfile,$category)=@ARGV;

#cortex_ac_v5_all_bypeak_enh_140723.bed
#phylop.vertebrate


my $filename;
if($bedfile=~/(.+).bed/) {
	$filename=$1;
}

#read in enh annotation #do calculation for CS16, CS23, F2 separately
my %clade2contents;
my %clade2name=qw(Human Human Ape Ape Primate Primate Amniota Amniota Gnathostomata MtAmn Mammalia Mammalia Tetrapoda MtAmn Vertebrata MtAmn Eutheria Eutheria Theria Theria MtAmn MtAmn NONE:Amniota NONE NONE:Tetrapoda NONE NONE:Gnathostomata NONE NONE:Vertebrata NONE NONE NONE);
my @names=qw(Human Ape Primate Eutheria Theria Mammalia Amniota MtAmn);

#annotation: sample->clade->peak
my %sample2peak;
my %enh2len;
open(IN,"/Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/MamCortex/ChIPseq/V5/Human/peaks_anno/cortex_ac_v5_all_annotation_combo_140723.txt") || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	next unless $array[1] eq "Enhancer";
	
	if($array[5] eq "embryonic_cortex") {
		foreach my $sample (split(";",$array[4])) {
			$sample2peak{$sample}{$clade2name{$array[6]}}{$array[0]}++;
			#$names{$clade2name{$array[6]}}++;
			if($array[0]=~/\w+\:(\d+)\-(\d+)/) {
				$enh2len{$array[0]}+=$2-$1;
			}
		}
	}
}
close IN;

#whole enhancer score
my $enhfile=$filename."_intersect_$category.txt";
my %enh2score;
open(IN,$enhfile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$enh2score{$array[3]}+=$array[8];
}
close IN;

#cce to len
my %enh2cce;
my $enh2cce_file=$filename."_intersect_ccemerged.bed";
open(IN,$enh2cce_file) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	my ($enh,$clade)=split(";",$array[3]);
	$enh2cce{$enh}{$clade2name{$clade}}+= $array[2]-$array[1];
}
close IN;

#cce to score
my $ccefile=$filename."_intersect_ccemerged_intersect_$category.txt";
my %cce2score;
open(IN,$ccefile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	my ($enh,$clade)=split(";",$array[3]);
	$cce2score{$enh}{$clade2name{$clade}}+=$array[8];
}
close IN;

#rpts to len
my %enh2rpt;
my $enh2rpt_file=$filename."_intersect_ccemerged_vs_rpts.bed";
open(IN,$enh2rpt_file) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	my ($enh,$clade,$rpt)=split(";",$array[3]);
	$enh2rpt{$enh}{$clade2name{$clade}}+= $array[2]-$array[1];
}
close IN;


#rpts to score
my $rptsfile=$filename."_intersect_ccemerged_vs_rpts_intersect_$category.txt";
my %rpts2score;
open(IN,$rptsfile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	my ($enh,$clade,$rpt)=split(";",$array[3]);
	$rpts2score{$enh}{$clade2name{$clade}}+=$array[8];
}
close IN;

print STDERR "Here!\n";

foreach my $sample (sort keys %sample2peak) {
	my $outfile=$filename."_$sample\_phastconsel_$category\_summary_forccenrpts.txt";
	open(OUT,">$outfile") || die $!;
	print OUT "Enhancer\tClade\t",join("\t",map {($_."_Len",$_."_LenOfElement")} qw(WholeEnhancer ConservedElement RptsInConservedElement)),"\n";
	foreach my $clade (@names) {
		foreach my $enh (sort keys %{$sample2peak{$sample}{$clade}}) {
			print OUT $enh,"\t",$clade,"\t";

			#whole enh
			print OUT $enh2len{$enh},"\t";
			if(defined $enh2score{$enh}) {
				print OUT $enh2score{$enh},"\t";
			}
			else {
				print OUT "NA\t";
			}
			
			#cce
			if(defined $enh2cce{$enh}{$clade}) {
				print OUT $enh2cce{$enh}{$clade},"\t";
			}
			else {
				print OUT "NA\t";
			}
			if(defined $cce2score{$enh}{$clade}) {
				print OUT $cce2score{$enh}{$clade},"\t";
			}
			else {
				print OUT "NA\t";
			}
			
			#rpt
			if(defined $enh2rpt{$enh}{$clade}) {
				print OUT $enh2rpt{$enh}{$clade},"\t";
			}
			else {
				print OUT "NA\t";
			}			
			if(defined $rpts2score{$enh}{$clade}) {
				print OUT $rpts2score{$enh}{$clade},"\n";
			}
			else {
				print OUT "NA\n";
			}
		}
	}
	close OUT;
}
	





