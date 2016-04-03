#!/usr/bin/perl -w
use strict;

my ($infile,$outfile)=@ARGV;

my @testlins=qw(Human Ape Primate Eutheria Theria Mammalia Amniota Tetrapoda Gnathostomata Vertebrata);
my @testspecs=qw(hg19	panTro2	ponAbe2	rheMac2	calJac1	mm9	rn4	cavPor3	oryCun2	bosTau4	equCab2	canFam2	loxAfr3	monDom5	ornAna1	galGal3	taeGut1	anoCar1	xenTro2	tetNig2	fr2	gasAcu1	oryLat2	danRer6	petMar1);
my @matchedlins=qw(Human Ape Ape Primate Primate Eutheria Eutheria Eutheria Eutheria Eutheria Eutheria Eutheria Eutheria Theria Mammalia Amniota Amniota Amniota Tetrapoda Gnathostomata Gnathostomata Gnathostomata Gnathostomata Gnathostomata Vertebrata);

#my %clade2col=(Human=>[1],Ape=>[2,3],Primate=>[4,5],Eutheria=>[6..13],Theria=>[14],Mammalia=>[15],Amniota=>[16..18],Tetrapoda=>[19],Gnathostomata=>[20..24],Vertebrata=>[25]);
    
#read processed age results
open(IN,$infile) || die $!;
my %peak2accumark;
my %peak2acculin;

while(<IN>) {
    tr/\r\n//d;
    my @array=split/\t/;
    next if $_=~/^Peak/;
    $peak2accumark{$array[0]}{$array[1]}=[@array[2..26]];
    $peak2acculin{$array[0]}{$array[1]}=[@array[27..36]];
}
close IN;

#print oUT results
open(OUT,">$outfile") || die $!;
print OUT "Peak\tCategory\t","MostDistantlyConsistentlyID\tMostDistantlyConsistentlySpec\tMostDistantlyConsistentlyClade\t","MostDistantlyID\tMostDistantlySpec\tMostDistantlyClade\t";
print OUT join("\t",@testspecs),"\t",join("\t",@testlins),"\n";
foreach my $peak (sort keys %peak2accumark) {
	foreach my $cate (sort keys %{$peak2acculin{$peak}}) {
			my $distantspec=decide_spec($peak2accumark{$peak}{$cate});
			my $distantspec_withinter=decide_spec_withinter($distantspec,$peak2acculin{$peak}{$cate});
			print OUT $peak,"\t",$cate,"\t";
			print OUT $distantspec_withinter,"\t",$distantspec_withinter==0?"NONE".":".$testspecs[$distantspec-1]:$testspecs[$distantspec-1],"\t",$distantspec_withinter==0?"NONE".":".$matchedlins[$distantspec-1]:$matchedlins[$distantspec_withinter-1],"\t";
			print OUT $distantspec,"\t",$testspecs[$distantspec-1],"\t",$matchedlins[$distantspec-1],"\t";
			print OUT join("\t",@{$peak2accumark{$peak}{$cate}}),"\t",join("\t",@{$peak2acculin{$peak}{$cate}}),"\n";
	}
}
close OUT;


#--------------------------------------------------------
sub decide_spec {
	my $specs=shift @_;
	my $spec=0;
	for(my $num=0;$num<@{$specs};$num++) {
		if($specs->[$num]>0) {
			$spec=$num;
		}
	}
	return $spec+1; #real num of most distantly aligned spec
}


sub decide_spec_withinter {
	my ($spec,$lins)=@_;
	my $specw=0;
	#start from the most conserved way, 
	#theria, with eutheria
	#mammal with eutheria
	#amniota with eutheria + mammalia
	#tetra with euteria  + mammalia + amiota
	#gnathos with eutheria + amiota + mam + tetra???
	#vete with eutheria + amiota + mam + tetra + gna
	my %lin2order=qw(Human 1 Ape 2 Primate 3 Eutheria 4 Theria 5 Mammalia 6 Amniota 7 Tetrapoda 8 Gnathostomata 9 Vertebrata 10);
	
	my $matchlin=$matchedlins[$spec-1];
	my $linno=0;
	if($lin2order{$matchlin}>5) {
		for(my $num=6;$num<=$lin2order{$matchlin};$num++) {
			if($lins->[$num-1]>0) {
				$linno++;
			}
		}
		if($linno==$lin2order{$matchlin}-5) {
			return $spec;
		}
		else {
			return 0;
		}
	}
	else {
		return $spec;
	}
}


