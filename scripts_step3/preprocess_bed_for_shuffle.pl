#!/usr/bin/perl -w
use strict;

my ($infile,$bedannotation, $background,$outfile) =@ARGV;

my %bed2coord;
#this has to be non-overlaping beds
open(IN,$infile) || die $!;
while(<IN>){
	tr/\r\n//d;
	my @array=split/\t/;
	$bed2coord{$array[3]}=[@array[0..2]];
}
close IN;

#bed to len/index
my %bed2anno;
my %lin2bed;
open(IN,$bedannotation) || die $!;
while(<IN>){
	tr/\r\n//d;
	my @array=split/\t/;
	next if $array[2] =~/NONE/;
	if(defined $bed2coord{$array[0]}) {
		$bed2anno{$array[0]}=[@array[2,6]];
		$lin2bed{$array[6]}{$bed2coord{$array[0]}[0]}{$array[2]}{$array[0]}=$bed2coord{$array[0]}[2]-$bed2coord{$array[0]}[1]; #lin, chr, reg, name => len
	}
}
close IN;

#background to len/index
my %cce2anno;
my %cce2bklen;
my %cce2background;
my %index2cce;
open(IN,$background) || die $!;
while(<IN>){
	tr/\r\n//d;
	my @array=split/\t/;
	$cce2anno{$array[3]}=[@array];
	$cce2background{$array[5]}{$array[0]}{$array[4]}{$array[3]}=$array[7]-$array[6]; #lin, chr,reg, name => len
	#$cce2bklen{$array[3]}{$array[5]}{$array[0]}{$array[4]}=$array[7]-$array[6];
	$index2cce{$array[5]}{$array[0]}{$array[4]}{$array[8]}=$array[3];
}
close IN;


#the output file should be sorted by length, with index???
open(OUT,">$outfile") || die $!;
my %bed2index;
my %bed2bind;
foreach my $lin (sort keys %lin2bed) {
	foreach my $chr (sort keys %{$lin2bed{$lin}}) {
		foreach my $reg (sort keys %{$lin2bed{$lin}{$chr}}) {
			#sort by len
			my $curind=1;
			my $bedindex=1;
			my $newreg="100ks";
			my $newcurind=1;
			my $newindex=1;			
			foreach my $bed (sort {$lin2bed{$lin}{$chr}{$reg}{$b} <=> $lin2bed{$lin}{$chr}{$reg}{$a}} keys %{$lin2bed{$lin}{$chr}{$reg}}) {
				#all bed files from a certain distance $reg, on the chrosome $chr, from a certain lineage $lin
				#$bed2index{$lin}{$chr}{$reg}{$bed}=$bedindex;
				
				while(defined $index2cce{$lin}{$chr}{$reg}{$curind}) {
					if($cce2background{$lin}{$chr}{$reg}{$index2cce{$lin}{$chr}{$reg}{$curind}}>=$lin2bed{$lin}{$chr}{$reg}{$bed}) {
						if(defined $index2cce{$lin}{$chr}{$reg}{$curind+1}) {
							$curind++;
						}
						else {
							last;
						}
					}
					else {
						$curind--; #go back to previous ID
						last;
					}
				}
				
				unless ($curind < $bedindex) {
					#$curind-1 is the last cce with longer shuffle region
					print OUT join("\t",$bed,$lin,$chr,$reg,$curind,$bedindex),"\n"; #enh, lin, chr,reg, backgroundindex, bedindex
				}
				else {	
					#not enough background reigon to shuffle
					if($reg eq "10k") {
						#use 100ks
						while(defined $index2cce{$lin}{$chr}{$newreg}{$newcurind}) {
							if($cce2background{$lin}{$chr}{$newreg}{$index2cce{$lin}{$chr}{$newreg}{$newcurind}}>=$lin2bed{$lin}{$chr}{$reg}{$bed}) {
								if(defined $index2cce{$lin}{$chr}{$newreg}{$newcurind+1}) {
									$newcurind++;
								}
								else {
									last;
								}
							}
							else {
								$newcurind--; #go back to previous ID
								last;
							}
						}
						unless ($newcurind < $newindex) {
							#$curind-1 is the last cce with longer shuffle region
							print OUT join("\t",$bed,$lin,$chr,$newreg,$newcurind,$newindex),"\n"; #enh, lin, chr,reg, backgroundindex, bedindex
						}
						else {
							print OUT join("\t",$bed,$lin,$chr,$newreg,$newcurind,$newindex),"\n"; #enh, lin, chr,reg, backgroundindex, bedindex
							print STDERR join("\t",$bed,$lin,$chr,$newreg,$newcurind,$newindex),"\n";
						}					
						$newindex++;
					}
					else {
						print OUT join("\t",$bed,$lin,$chr,$reg,$curind,$bedindex),"\n"; #enh, lin, chr,reg, backgroundindex, bedindex
						print STDERR join("\t",$bed,$lin,$chr,$reg,$curind,$bedindex),"\n"; #usually should not happen
					}
				}
				$bedindex++;
			}
		}
	}
}
close OUT;
				
				
					
			



