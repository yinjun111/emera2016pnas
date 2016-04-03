#!/usr/bin/perl -w
use strict;

#read in enh annotation #do calculation for CS16, CS23, F2 separately
my %clade2contents;
my %clade2name=qw(Human Human Ape Ape Primate Primate Amniota Amniota Gnathostomata MtAmn Mammalia Mammalia Tetrapoda MtAmn Vertebrata MtAmn Eutheria Eutheria Theria Theria NONE:Amniota NONE NONE:Tetrapoda NONE NONE:Gnathostomata NONE NONE:Vertebrata NONE NONE NONE MtAmn MtAmn);
my @names=qw(Human Ape Primate Eutheria Theria Mammalia Amniota MtAmn);

#annotation: sample->clade->peak
my %sample2peak;
open(IN,"/Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/MamCortex/ChIPseq/V5/Human/peaks_anno/cortex_ac_v5_all_annotation_combo_140723.txt") || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	next unless $array[1] eq "Enhancer";
	
	if($array[5] eq "embryonic_cortex") {
		foreach my $sample (split(";",$array[4])) {
			$sample2peak{$sample}{$clade2name{$array[6]}}{$array[0]}++;
			#$names{$clade2name{$array[6]}}++;
		}
	}
}
close IN;


#intersection
my %peak2len;
open(IN,"cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_vs_rpts.txt") || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	my ($peak,$clade)=split(";",$array[3]);
	$peak2len{$peak}{$clade2name{$clade}}+=$array[8];
}
close IN;

#result
open(OUT,">cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_vs_rpts_lensummary.txt") || die $!;
print OUT "Stage\t",join("\t",@names),"\n";
foreach my $sample (sort keys %sample2peak) {
	my @nums;
	foreach my $name (@names) {
		if(defined $sample2peak{$sample}{$name}) {
			my $len=0;
			foreach my $peak (sort keys %{$sample2peak{$sample}{$name}}) {
				if(defined $peak2len{$peak}{$name}) {
					$len+=$peak2len{$peak}{$name};
				}
			}
			push @nums,$len;
		}
		else {
			push @nums,0;
		}
	}
	print OUT $sample,"\t",join("\t",@nums),"\n";
}
close OUT;
