#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);

#compute give results line by line
#summarize all results from WIN LEN
#result is the accumulated number of windows > threshold

my @testlins=qw(Human Ape Primate Eutheria Theria Mammalia Amniota Tetrapoda Gnathostomata Vertebrata);
my @testspecs=qw(hg19	panTro2	ponAbe2	rheMac2	calJac1	mm9	rn4	cavPor3	oryCun2	bosTau4	equCab2	canFam2	loxAfr3	monDom5	ornAna1	galGal3	taeGut1	anoCar1	xenTro2	tetNig2	fr2	gasAcu1	oryLat2	danRer6	petMar1);

my @selcols=3..27;

my ($infile,$outfile,$perc)=@ARGV;
#raw age, processed age, 0.5

open(IN,$infile) || die $!;
my %peak2accumark;
my %peak2acculin;

while(<IN>) {
    tr/\r\n//d;
    my @array=split/\t/;

    if ($_=~/^File/ ) {
    	next;
    }
    #name of the enhancer
	my $peak;
	if($array[0]=~/([^\/]+)\.fas$/) {
		$peak=$1;
	}

    #error. shouldn't happen
    if($array[$selcols[0]]==0) {
        print STDERR $_,"\n";
        next;
    }
    
    #decide status of each row
    my @marks=mark_div(1,$perc,@array[@selcols]);
    #decide status of the lins
    my @linmarks=decide_lin(@marks);
	
	my $cate;
	if($array[2]=~/^Whole/) {
		$cate=$array[1]."_Whole";
	}
	else {
		$cate=$array[1]."_Win";
	}
	
	if(defined $peak2accumark{$peak}{$cate}) {
		$peak2accumark{$peak}{$cate}=add_mark(\@{$peak2accumark{$peak}{$cate}},\@marks);
	}
	else {
		$peak2accumark{$peak}{$cate}=[@marks];
	}

	if(defined $peak2acculin{$peak}{$cate}) {
		$peak2acculin{$peak}{$cate}=add_mark(\@{$peak2acculin{$peak}{$cate}},\@linmarks);
	}
	else {
		$peak2acculin{$peak}{$cate}=[@linmarks];
	}
}

open(OUT,">$outfile") || die $!;
print OUT "Peak\tCategory\t",join("\t",@testspecs),"\t",join("\t",@testlins),"\n";
foreach my $peak (sort keys %peak2accumark) {
	foreach my $cate (sort keys %{$peak2accumark{$peak}}) {
		print OUT $peak,"\t",$cate,"\t";
		print OUT join("\t",@{$peak2accumark{$peak}{$cate}}),"\t";
		print OUT join("\t",@{$peak2acculin{$peak}{$cate}}),"\n";
	}
}
close OUT;


#--------------

sub mark_div {
	#give mark on each species based on peak alignment
    my ($ref,$perc,@lens)=@_;
    my @marks;
    for(my $num=0;$num<@lens;$num++) {
	if($lens[$ref-1]>0 && $lens[$num]>0 && $lens[$num]/$lens[$ref-1]>=$perc) {
	    $marks[$num]=1;
	}
	else {
	    $marks[$num]=0;
	}
    }
    return @marks;
}

sub decide_lin {
	#read mark from each spec, and assign mark to lineage
    my @marks=@_;
    my @linmarks;
	my %clade2col=(Human=>[1],Ape=>[2,3],Primate=>[4,5],Eutheria=>[6..13],Theria=>[14],Mammalia=>[15],Amniota=>[16..18],Tetrapoda=>[19],Gnathostomata=>[20..24],Vertebrata=>[25]);
    
    foreach my $lin (@testlins) {
    	my @specmarks=@marks[map {$_-1} @{$clade2col{$lin}}];
    	push @linmarks,sum(@specmarks);
    }
    return @linmarks; #number of specs matched in each lineage
}

sub add_mark {
	my ($array1,$array2)=@_;
	my @results;
	for(my $num=0;$num<@{$array1};$num++) {
		$results[$num]=$array1->[$num]+$array2->[$num];
	}
	return [@results];
}

#----------------------------------------

