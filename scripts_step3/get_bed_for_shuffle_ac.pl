#!/usr/bin/perl -w
use strict;

my %cate2value=("Human",1,"Ape",2,"Primate",3,"Eutheria",4,"Theria",5,"Mammalia",6,"Amniota",7,"Tetrapoda",8,"Gnathostomata",9,"Vertebrata",10,"NONE",0,"UNDEF",-1);

open(IN,"/Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/Repeat/ChIPseq/V5/Human/peaks_anno/cortex_ac_v5_all_annotation_combo_140723.txt") || die $!;
open(OUT,">cortex_ac_v5_all_bypeak_enh_140723.bed") || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	#Enhancer, with match to mouse, and more conserved than Theria
	if($array[1] eq "Enhancer" && $array[5] ne "None" && $array[5] ne "Nomatch") {
		if($cate2value{$array[6]}>=4) {
			#mammal and beyond
			if($array[0]=~/(\w+)\:(\d+)\-(\d+)/) {
				print OUT $1,"\t",$2,"\t",$3,"\t",$array[0],"\n";
			}
		}
	}
}
close IN;
close OUT;
