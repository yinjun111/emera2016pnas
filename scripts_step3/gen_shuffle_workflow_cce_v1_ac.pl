#!/usr/bin/perl -w
use strict;

#Get a general purpose script for all sorts of shuffling
#A uniformed shuffling for all analyses

#------------------------------
#generate shuffle, and count it
#general info
my $shuffletime=20000;
my $bedfile="cortex_ac_v5_all_bypeak_enh_140723.bed"; #name for the original shuffle bedfile
my $annofile="cortex_ac_v5_all_bypeak_enh_140723_anno_forshuffle.txt";
my $bkfile="gencode_v10_cce_all_shuffle_background_forShuffle_140606.txt";

my $filename;
if($bedfile=~/([^\/]+)\.bed/) {
	$filename=$1;
}

#more general parameters
my $shuffledir="/home/jy344/Scratch/Projects/Repeats/shuffled/shuffle_cce/ac";

#use automatic info for the other files
#pbs files
my $shufflepbs="$filename\_shuffle$shuffletime.pbs";

#------------------------------
#!!! Now it works
shuffle_pbs($bedfile,$filename,$shuffledir,$shuffletime,$shufflepbs);

#Functions for the shuffle workflow
#------------------------------
sub shuffle_bed {
	my ($infile,$outfile)=@_; #only need to consider in and out
	#whatever shuffling method to use
	#Global: $annofile $bkfile #other parameters depending on the methods
	return "perl /home/jy344/Programs/Scripts/Rmsk/shuffle_bed_intersect_cce_140102.pl $infile $annofile $bkfile $outfile";
	#return "/home/jl56/TOOLS/BEDTools-Version-2.10.0/bin/shuffleBed -i $infile -g $chrom -excl $exclfile > $outfile"; #an alternative for CorDev
}


#------------------------------
#uniform shuffle, count and sum workflow

sub shuffle_pbs {
	my ($bedfile,$filename,$workdir,$shuffletime,$outfile)=@_;
	#Global: none
	#Command for shuffled bed
	#------------------------------
	open(OUT,">$outfile") || die $!;
	foreach my $num (1..$shuffletime) {
		print OUT "cd $workdir;";
		#shuffle & intersect
		print OUT shuffle_bed($bedfile,"$filename\_s$num.bed");
		print OUT "\n";
	}
	close OUT;
}
