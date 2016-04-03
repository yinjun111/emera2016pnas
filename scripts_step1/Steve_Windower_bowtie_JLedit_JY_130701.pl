#!/usr/bin/perl 

#perl ~/Scripts/Windower_bowtie.pl e11.5_limb e11.5_limb_input ~/annotations/mouse/mm9/chrom.sizes .001 -v
use warnings;
use strict;
use lib "/home/jy344/Programs/perl/lib/lib/perl5/x86_64-linux-thread-multi/"; #JY
use lib "/home/jy344/Programs/perl/lib/lib/perl5/"; #JY
use Math::CDF;
use Parallel::ForkManager;#JY

my %chrom_sizes= ();
#my $ChIP_wig= $ARGV[0]."ChIP.wig";
#my $input_wig= $ARGV[0]."input.wig";
my $pval= $ARGV[3];
my $out_bed= $ARGV[0]."EnrichedRegions_temp.bed";
my $out_verb= $ARGV[0]."Data-riffic".$pval.".txt";

#JY
unless(-r $ARGV[0]."tmp") {
	mkdir($ARGV[0]."tmp");
}

my $cpu=$ARGV[5];

#JY
unless(defined $cpu) {
	$cpu=8;
}

my $pm; #object for parallel
unless($cpu == 1) {
	$pm=Parallel::ForkManager->new($cpu,$ARGV[0]."tmp"); #JY #Edit number of threads you want to use
}


#deal with wig, verbose arguments
#my $wiggin = 'nope';
#if (defined $ARGV[5]){
#	$wiggin = $ARGV[5];
#}

my $verbosin = 'nope';
if (defined $ARGV[4]){
	$verbosin = $ARGV[4];
}

#if ("-wig" eq $wiggin){
#	open (CHIP,">$ChIP_wig");
#	open (INPUT,">$input_wig");
#	print "This time, we're wiggin' out. \n";
#}

if ("-v" eq $verbosin){
	open (VERB,">$out_verb");
	print "This time, you're getting so much extra info. \n";
	print VERB "Chrom\twindow_start\tchipTags\tinputTags\tpval \n";
}




my $Chip_input=$ARGV[0];
my $Input_input=$ARGV[1];
my $Chrom_input=$ARGV[2];

#check to make sure chom.sizes files exists
if (!open(INPUT3,$Chrom_input)){
	print STDERR "Can't open file: ".$Chrom_input." Make sure its there.\n";
	exit;
   	}
#get chromosomes and sizes   	
while(my $tmp1=<INPUT3>){
	chomp ($tmp1);
	my ($temp_chrom,$temp_End)=split(/\s+/,$tmp1);
	$chrom_sizes{$temp_chrom}=$temp_End;
}

my @chrs1= sort keys %chrom_sizes; 
my @chrs= grep(!/rand/, @chrs1);

#figure out genome size
my $TotGenSize=0;
foreach my $chr (@chrs){
	$TotGenSize= $TotGenSize + $chrom_sizes{$chr};
}

print "The calculated genome size is ", $TotGenSize,"\n";

#check to make sure sam files exists
  	if (!open(INPUT1,$Chip_input.".aligned")){
		print STDERR "Can't open file: ".$Chip_input." Make sure its there, also take off .sam.\n";
		exit;
    	}  
my $test=<INPUT1>; #JLedit
chomp ($test); #JLedit
my @test_split=split(/\s+/,$test); #JLedit
my $Tag_len=length ($test_split[5]); #JLedit
print ("ChIP tag length: ",$Tag_len,"\n"); #JLedit

	if (!open(INPUT2,$Input_input.".aligned")){
		print STDERR "Can't open file: ".$Input_input." Make sure its there, also take off .sam.\n";
		exit;
    	}   
$test=<INPUT2>; #JLedit                                                                                          
chomp ($test); #JLedit                                                                                              
@test_split=split(/\s+/,$test); #JLedit                                                                          
my $IN_len=length ($test_split[5]); #JLedit                                                                        
print ("Input tag length: ",$IN_len,"\n"); #JLedit  

print ("Preppin' and cleanin' your sam files. \n");

#system("samtools view -Sq 1 ". $Chip_input. ".sam > ". $Chip_input. "_reliable.sam");
#system("samtools view -Sq 1 ". $Input_input. ".sam > ". $Input_input. "_reliable.sam");
system("cut -f 2,3,4 ". $Chip_input. ".aligned > ". $Chip_input. "_cut.sam");
system("cut -f 2,3,4 ". $Input_input. ".aligned > ". $Input_input. "_cut.sam");
#system("rm ". $Chip_input. "_reliable.sam");
#system("rm ". $Input_input. "_reliable.sam");
close INPUT1;
close INPUT2;

print ("Zawesome, begining magic. \n");

#iterate over each chromosome	
foreach my $chr (@chrs){
	unless($cpu == 1) {
		$pm->start and next;#JY
	}
	print $chr, "\n";

	open (BED, ">${out_bed}_$chr"); #JY

# my %tag_orientations=(); JLedit
# my %IN_tag_orientations=(); JLedit
 
 my %Tag_middle=();
 my %IN_middle=();
 
 my %EnrichedWindows =();
 
	#check to make sure sam files exists
  	if (!open(INPUT1,$Chip_input."_cut.sam")){
		print STDERR "Can't open file: ".$Chip_input." Make sure its there, also take off .sam.\n";
		exit;
    	}   
	if (!open(INPUT2,$Input_input."_cut.sam")){
		print STDERR "Can't open file: ".$Input_input." Make sure its there, also take off .sam.\n";
		exit;
    	}   

	#Create 25bp markers and Window/Bins
	 print "Binnin' \n";
 	 my %Bin_Ranges=();
 	 my %IN_Bin_Ranges=();
 	 my %Wind_Ranges=();
 	 my %IN_Wind_Ranges=();
 	 
 	 for (my $i=1; $i <= $chrom_sizes{$chr};){
 		 $Bin_Ranges{$i}=0;
 		 $IN_Bin_Ranges{$i}=0;
 		 $Wind_Ranges{$i}=0;
 		 $IN_Wind_Ranges{$i}=0;
 		 $i=$i+25;
 	 }
 
	#sort bin keys and produce array with names of bins.  
	my @bin_keys = sort {$a<=>$b} keys %Bin_Ranges; 
	print "Done Binnin', sorry that took so long, the next one will be faster, I promise.\n";


#input Sam alignments by chromosome
    my $tag_count=0;
    while(my $tmp1=<INPUT1>){
	chomp ($tmp1);
	my ($flag,$chrm,$T_start)=split(/\s+/,$tmp1);
	#pull out lines that match chromosome, make hash with tag_start as key and orientation flag as value
	if ($chrm eq $chr){
	    if ($flag eq "+") {$Tag_middle{($T_start+150).$flag}=$flag;} #JLedit 
	    else {$Tag_middle{($T_start-150+$Tag_len-1).$flag}=$flag;} #JLedit 
	}
	$tag_count=scalar keys %Tag_middle; #JLedit
    }
       
    print "ChIP Tag Countin', ", $tag_count, " tags.\n";   
    
    my $IN_count=0;
    while(my $tmp2=<INPUT2>){
	chomp ($tmp2);
	my ($flag,$chrm,$T_start)=split(/\s+/,$tmp2);
	#pull out lines that match chromosome, make hash with tag_start as key and orientation flag as value
	if ($chrm eq $chr){
            if ($flag eq "+") {$IN_middle{($T_start+150).$flag}=$flag;} #JLedit                                    
	    else {$IN_middle{($T_start-150+$IN_len-1).$flag}=$flag;} #JLedit    
	}
	$IN_count=scalar keys %IN_middle; #JLedit
    }    
    print "Input Tag Countin'n', ", $IN_count," tags.\n";  


	#create histogram in 25bp blocks and create hashes of inferred tag centers as key and orig tstart as value
	
	# for chip tags
#	my @Tag_keys = sort {$a<=>$b} keys %tag_orientations;
#  	foreach my $Tag_key (@Tag_keys){
#		if ($tag_orientations{$Tag_key} eq '+'){
#			$Tag_middle{$Tag_key+150}=$Tag_key;
#			if ("-wig" eq $wiggin){
#				for (my $i=$Tag_key; $i<=($Tag_key+300); $i++){
#					if (exists $Bin_Ranges{$i}){
#						$Bin_Ranges{$i}= $Bin_Ranges{$i} + 1;					
#					}
#				}
#			}		
#		}
#		if ($tag_orientations{$Tag_key} eq '-'){
#			$Tag_middle{$Tag_key-150}=$Tag_key;
#			if ("-wig" eq $wiggin){
#				for (my $i=($Tag_key-300); $i<=$Tag_key; $i++){
#					if (exists $Bin_Ranges{$i}){
#						$Bin_Ranges{$i}= $Bin_Ranges{$i} + 1;	
#					}
#				}
#			}
#		}	
#	}
	
#maybe clear tag hashes here to save memory
	
	# for input tags
#	my @IN_Tag_keys = sort {$a<=>$b} keys %IN_tag_orientations;
#  	foreach my $IN_Tag_key (@IN_Tag_keys){
#		if ($IN_tag_orientations{$IN_Tag_key} eq '+'){
#			$IN_middle{$IN_Tag_key+150}=$IN_Tag_key;
#			if ("-wig" eq $wiggin){
#				for (my $i=($IN_Tag_key-300); $i<=$IN_Tag_key; $i++){
#					if (exists $Bin_Ranges{$i}){
#						$IN_Bin_Ranges{$i}= $IN_Bin_Ranges{$i} + 1;
#					}
#				}
#			}
#		}
#		if ($IN_tag_orientations{$IN_Tag_key} eq '-'){
#			$IN_middle{$IN_Tag_key-150}=$IN_Tag_key;
#			if ("-wig" eq $wiggin){
#				for (my $i=$IN_Tag_key; $i<=($IN_Tag_key-300); $i++){
#					if (exists $Bin_Ranges{$i}){
#						$IN_Bin_Ranges{$i}= $IN_Bin_Ranges{$i} + 1;	
#					}
#				}
#			}
#		}
#	}

	# Do crazy statistics
	print "Doin' Stats (Head Hurtin').\n";
    my $lamda_tags=  (500 * $tag_count / $chrom_sizes{$chr}/0.85);  #JLedit
    my $IN_lamda_tags=  (500 * $IN_count / $chrom_sizes{$chr}/0.85); #JLedit
	
	
	#Generate Bin Counts
	foreach my $bin_key (@bin_keys){
		if ($bin_key % 10000000 ==1){
			my $percent = $bin_key / $chrom_sizes{$chr} *100;
			print $percent, " % done hurtin' \n";
		}
		for (my $i=$bin_key; $i<$bin_key+500;$i++){ #JLedit
			
			if (exists $Tag_middle{$i."+"}){ #JLedit
				$Wind_Ranges{$bin_key}=$Wind_Ranges{$bin_key}+1;
			}
			if (exists $Tag_middle{$i."-"}){ #JLedit
			        $Wind_Ranges{$bin_key}=$Wind_Ranges{$bin_key}+1;
			}
			if (exists $IN_middle{$i."+"}){ #JLedit
				$IN_Wind_Ranges{$bin_key}=$IN_Wind_Ranges{$bin_key}+1;
			}
                        if (exists $IN_middle{$i."-"}){ #JLedit                                                 
			    $IN_Wind_Ranges{$bin_key}=$IN_Wind_Ranges{$bin_key}+1;
                        }
		}
		
		if ($Wind_Ranges{$bin_key} >= 1){
		
				my $epsilon= $IN_Wind_Ranges{$bin_key} / $IN_lamda_tags;
				my $max_choice = 1;
				if ($epsilon > $max_choice){
					$max_choice = $epsilon;
					#print "i make choices, ",$max_choice," \n";
				}
				my $lamda = $max_choice * $lamda_tags;
				my $pVal=(1 - &Math::CDF::ppois($Wind_Ranges{$bin_key}, $lamda));
				if ($pVal <= $pval){
					$EnrichedWindows{$bin_key}=$pVal;
					if ("-v" eq $verbosin){
						print VERB $chr, "\t", $bin_key, "\t", $Wind_Ranges{$bin_key}," \t", $IN_Wind_Ranges{$bin_key},  "\t",$pVal, "\n";
						#print "added! \n";
					}
				}
				#print $bin_key, "\t", $epsilon, "\t", $max_choice, "\t", $lamda, "\t",$lamda_tags, "\t", $Wind_Ranges{$bin_key}," \t", $IN_Wind_Ranges{$bin_key},  "\t",$pVal, "\t", $tag_count, "\t",$IN_count, "\n",;			
		}
	}
	
	print "Makin' Yo' .Bed.\n";
	my @rich_keys = sort {$a<=>$b} keys %EnrichedWindows;
	print "Regions in $chr :",scalar(@rich_keys), "\n"; #JY
  	foreach my $rich_key (@rich_keys){
  		my $end_site=$rich_key +500; 
  		print BED $chr, "\t", $rich_key-1, "\t",  $end_site-1, "\t", $EnrichedWindows{$rich_key}, "\n"; #JLedit
  	}
  	
#  	if ("-wig" eq $wiggin){
  	
#    	print "Fixin' Yo' .Wigs.\n";
    	#make wig file for chip and input
#    	print CHIP "fixedStep  chrom=", $chr, " start=1  step=25  span=25 \n";
#    	foreach my $bin_key (@bin_keys){
#    		my $scratch = $Bin_Ranges{$bin_key};
#    		my $normalized_binCount = $scratch * 10000000 / $tag_count;
#    		print CHIP $normalized_binCount, "\n";
#    	}
#    	print CHIP "\n";
  	
#  		print INPUT "fixedStep  chrom=", $chr, " start=1  step=25  span=25 \n";
#  		foreach my $bin_key (@bin_keys){
#  	  		my $scratch = $IN_Bin_Ranges{$bin_key};
#  	 		my $IN_normalized_binCount = $scratch * 10000000 / $IN_count;
#  	 		print INPUT $IN_normalized_binCount, "\n";
#  		}
#  	    print INPUT "\n";
    

    
#	}
	unless($cpu == 1) {
		$pm->finish; #JY
	}
}

unless($cpu == 1) {
	$pm->wait_all_children;#JY
}

print "ALMOST DONE(thank Darwin)!\n";
system("cat  $out_bed\_chr* | /home/jl56/TOOLS/BEDTools-Version-2.10.0/bin/mergeBed -d 1000 -i stdin > ". $ARGV[0]. "EnrichedRegions".$pval.".bed"); #JY
system("rm $out_bed\_chr*");
system("rm ". $Chip_input. "_cut.sam");
system("rm ". $Input_input. "_cut.sam");
close INPUT1;
close INPUT2;

#if ("-wig" eq $wiggin){	
#	close CHIP;
#	close INPUT; 
#}
if ("-v" eq $verbosin){
	close VERB;
}
print "Done-zo \n";



