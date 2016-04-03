#!/usr/bin/perl -w

#two ways to use the script. Either use it in the Script, if with default setting. Or copy it to the working folder, and customize your parameters

#parameters to choose
windowsize=10000
shortwindow=5000
mafwin=$windowsize #10000 or 1000000

ageana=yes
phastana=no

#create working directories
for workdir in beds mafs fas selfas age_results rphast_results rphast_beds rphast_trees
do
  if ! [ -e $workdir ]
      then
      mkdir $workdir
  fi
done


#read in bed file
infile=$1
pbsfile=${infile/.bed/_submit.pbs}

if [ -e $pbsfile ]
    then
    rm $pbsfile
fi


#anno element to window
perl /home/jy344/Programs/Scripts/Age/bin/edit_bed_for_window.pl $infile ${infile/.bed/_forwindow.bed}

perl /home/jy344/Programs/Scripts/Age/bin/get_region_window_131121.pl ${infile/.bed/_forwindow.bed} ${infile/.bed/_windowanno.txt} $windowsize $shortwindow

#split bed into seperate files
#perl /home/jy344/Programs/Scripts/Age/bin/split_bed_121016.pl $infile beds
directory=`pwd`

for window in `cut -f 2 ${infile/.bed/_windowanno.txt} | sort | uniq`;
	do
		echo "cd $directory; for bed in \`grep $window ${infile/.bed/_windowanno.txt} | cut -f 1 | sort | uniq\`;do sh /home/jy344/Programs/Scripts/Age/bin/cal_age_for_bed_131121.sh $infile $bed $window $ageana $phastana $mafwin;done;" >> $pbsfile
done;
