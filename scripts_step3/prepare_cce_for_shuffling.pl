#!/usr/bin/perl -w
use strict;

#my %cate2value=("NONE",0,"Human",1,"Ape",2,"Primate",3,"Eutheria",4,"Theria",5,"Mammalia",6,"Amniota",7,"Tetrapoda",8,"Gnathostomata",9,"Vertebrata",10);

my @cates=("Human","Ape","Primate","Eutheria","Theria","Mammalia","NONE:Amniota","Amniota","NONE:Tetrapoda","Tetrapoda","NONE:Gnathostomata","Gnathostomata","NONE:Vertebrata","Vertebrata");

my @regionfiles=glob("regions_merged/*.bed");

open(OUT,">gen_region_forcce.sh") || die $!;
for(my $num1=0;$num1<@cates;$num1++) {
	next if $cates[$num1]=~/^NONE/;
	my $infile="../merged/hg19_nh.chrom.sizes_tobins100_50_v65_".$cates[$num1]."_merged.bed";
	my @removed; #remove all regions more conserved than the current one
	for(my $num2=$num1+1;$num2<@cates;$num2++) {
		#remove all more conserved elements from the regions
		push @removed, "../merged/hg19_nh.chrom.sizes_tobins100_50_v65_".$cates[$num2]."_merged.bed";
	}
	
	foreach my $regionfile (@regionfiles) {
		my $outfile;
		if($regionfile=~/([^\/]+)\.bed/) {
			my $outfile1=$1."_".$cates[$num1]."_region.bed";
			my $outfile2=$1."_".$cates[$num1]."_intersect.txt";
			if(@removed) {
				print OUT "cat ",join(" ", @removed),"| subtractBed -a $regionfile -b stdin | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$1\":\"\$2\"-\"\$3}' > intersect/$outfile1;";
			}
			else {
				#for Vertebrata, it doesn't need to remove anything
				print OUT "cat $regionfile | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$1\":\"\$2\"-\"\$3}' > intersect/$outfile1;";
			}
			print OUT "intersectBed -a $infile -b intersect/$outfile1 -wa -wb > intersect/$outfile2\n";
		}
	}
}
close OUT;
