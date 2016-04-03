#!/bin/sh

infile1=$1
infile2=$2

#param
spec1="hg19"
spec2="mm9"
liftOverfile1=$3
liftOverfile2=$4

#
resultfile=$5

#threshold=0.5
threshold=$6

#liftOverfile1="/home/jy344/Genomes/alignments/hg19ToPanTro3.over.chain.gz"
#liftOverfile2="/home/jy344/Genomes/alignments/panTro3ToHg19.over.chain.gz"

#LiftOver human to mouse
#liftOver
/home/jy344/Programs/bin/liftOver -minMatch=0.00001 -multiple $infile1 $liftOverfile1 ${infile1/.bed/_to$spec2.bed} ${infile1/.bed/_failto$spec2.bed}
#Unique liftOver
perl /home/jy344/Programs/Scripts/Bed/liftOver_uniq.pl ${infile1/.bed/_to$spec2.bed} ${infile1/.bed/_to$spec2\_uniq.bed}
#no-random
grep -v "_" ${infile1/.bed/_to$spec2\_uniq.bed} > ${infile1/.bed/_to$spec2\_uniq_norandom.bed}
#intersect
if [ "$threshold" -lt 0 ]; then
    /home/jl56/TOOLS/BEDTools-Version-2.10.0/bin/intersectBed -a ${infile1/.bed/_to$spec2\_uniq_norandom.bed} -b $infile2 -wa -wb -f $threshold -r | cut -f 4,8 | sort | uniq > ${infile1/.bed/_to$spec2\_uniq_norandom_intersect$spec2.bed}
else 
    echo "here"
    /home/jl56/TOOLS/BEDTools-Version-2.10.0/bin/intersectBed -a ${infile1/.bed/_to$spec2\_uniq_norandom.bed} -b $infile2 -wa -wb | cut -f 4,8 | sort | uniq > ${infile1/.bed/_to$spec2\_uniq_norandom_intersect$spec2.bed}
fi



#LiftOver mouse to human
#liftOver
/home/jy344/Programs/bin/liftOver -minMatch=0.00001 -multiple $infile2 $liftOverfile2 ${infile2/.bed/_to$spec1.bed} ${infile2/.bed/_failto$spec1.bed}
#Unique liftOver
perl /home/jy344/Programs/Scripts/Bed/liftOver_uniq.pl ${infile2/.bed/_to$spec1.bed} ${infile2/.bed/_to$spec1\_uniq.bed}
#no-random
grep -v "_" ${infile2/.bed/_to$spec1\_uniq.bed} > ${infile2/.bed/_to$spec1\_uniq_norandom.bed}
#intersect
if [ "$threshold" -lt 0 ]; then
    /home/jl56/TOOLS/BEDTools-Version-2.10.0/bin/intersectBed -a ${infile2/.bed/_to$spec1\_uniq_norandom.bed} -b $infile1 -wa -wb -f $threshold -r | cut -f 4,8 | sort | uniq | awk '{print $2"\t"$1}' > ${infile2/.bed/_to$spec1\_uniq_norandom_intersect$spec1.bed}
else
    echo "here"
    /home/jl56/TOOLS/BEDTools-Version-2.10.0/bin/intersectBed -a ${infile2/.bed/_to$spec1\_uniq_norandom.bed} -b $infile1 -wa -wb | cut -f 4,8 | sort | uniq | awk '{print $2"\t"$1}' > ${infile2/.bed/_to$spec1\_uniq_norandom_intersect$spec1.bed}
fi


#find double match
cat ${infile1/.bed/_to$spec2\_uniq_norandom_intersect$spec2.bed}  ${infile2/.bed/_to$spec1\_uniq_norandom_intersect$spec1.bed} | sort | uniq -d > $resultfile
