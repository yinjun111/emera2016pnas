#!/usr/bin/perl -w
use strict;

my ($bedfile,$annofile,$ccefile,$outfile,$oloption) = @ARGV;
#bedfile in a certain clade from a defined region (distance to gene), intersect with cce from the same clade in the same defined region

unless(defined $oloption) {
	$oloption=0; #don't use cce again
}

#4 col bed file
my %bed2coord;
open(IN,$bedfile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$bed2coord{$array[3]}=[@array[0..2]];
}
close IN;

#bed file anotation
#bed, chr,clade, genomeregion,cceindex
my %bed2anno;
open(IN,$annofile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$bed2anno{$array[0]}=[@array[1..4]]; #lin, chr, reg, index
}
close IN;

#ccefile
#6 col bed-like file
#chr, start, end, cce, genomeregion,clade, estart, eend,cceindex
my %cce2info;
my %cate2cce;
open(IN,$ccefile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$cce2info{$array[3]}{$array[5]}{$array[0]}{$array[4]}=[@array]; #cce,lin, chr,reg
	$cate2cce{$array[5]}{$array[0]}{$array[4]}{$array[8]}=$array[3]; #lin, chr,reg,index,cce
}
close IN;

#shuffle
my %used_cce;
my $sel_cce;
open(OUT,">$outfile") || die $!;
foreach my $bed (sort keys %bed2coord) {
	$sel_cce=get_rand_cce(@{$bed2anno{$bed}},$oloption,\%used_cce);
	unless($sel_cce eq "None") {
		$used_cce{$sel_cce}++;
		my ($rand_bed_start,$rand_bed_end)=rand_by_cce($bed2coord{$bed},[@{$cce2info{$sel_cce}{$bed2anno{$bed}[0]}{$bed2anno{$bed}[1]}{$bed2anno{$bed}[2]}}[0,1,2,6,7]]);
		print OUT join("\t",$bed2coord{$bed}[0],$rand_bed_start,$rand_bed_end,$bed),"\n";
	}
	else {
		print STDERR "$bed no space for shuffling! Use dup region to shuffle!\n";
		$sel_cce=get_rand_cce(@{$bed2anno{$bed}},1);
		my ($rand_bed_start,$rand_bed_end)=rand_by_cce($bed2coord{$bed},[@{$cce2info{$sel_cce}{$bed2anno{$bed}[0]}{$bed2anno{$bed}[1]}{$bed2anno{$bed}[2]}}[0,1,2,6,7]]);
		print OUT join("\t",$bed2coord{$bed}[0],$rand_bed_start,$rand_bed_end,$bed),"\n";
	}
}
close OUT;


sub get_rand_cce {
	my ($lin,$chr,$reg,$cceindex,$ol,$used_cce)=@_;
	if ($ol==1) {
		#dup cce allowed
		return $cate2cce{$lin}{$chr}{$reg}{int(rand($cceindex)+1)};
	}
	else {
		#dup not allowed
		#get the not used list

		my @cand_index;
		foreach my $index (1..$cceindex) {
			if(defined $cate2cce{$lin}{$chr}{$reg}{$index}) {
				unless(defined $used_cce->{$cate2cce{$lin}{$chr}{$reg}{$index}}) {
					push @cand_index,$index;
				}
			}
			else {
				print STDERR join("\t",$lin,$chr,$reg,$index),"\n";
			}
		}
			
		unless(@cand_index==0) {
			#should not happen, no more cce to use
			return $cate2cce{$lin}{$chr}{$reg}{$cand_index[int(rand(@cand_index))]};
		}
		else {
			return "None";
		}
	}
}

sub rand_by_cce {
	my ($bedcoord,$ccecoord)=@_;
	#bedcoord: chr, start, end
	#ccecoord: chr, cstart, cend, bstart, bend
	
	#first determine left and right most coord
	#start, end are in BED format
	my ($pstart,$pend); #possible range for the start of the rand bed, both 0 based
	my $extralen=$bedcoord->[2] - $bedcoord->[1] -100; #extended len of peak to cce, at least overlap by 100bp
	if($ccecoord->[1] - $ccecoord->[3] >=$extralen) {
		$pstart=$ccecoord->[1] - $extralen;
	}
	else {
		$pstart=$ccecoord->[3];
	}
	if($ccecoord->[4] - $ccecoord->[2]>=$extralen) {
		$pend=$ccecoord->[2] -100;
	}
	else {
		$pend=$ccecoord->[4]-($bedcoord->[2] - $bedcoord->[1]);
	}
	my $rand_start=$pstart+int(rand($pend-$pstart+1));
	
	return ($rand_start,$rand_start+$bedcoord->[2] - $bedcoord->[1]); #in BED format
}		


