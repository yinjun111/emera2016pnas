#!/bin/sh

bedfile=$1
result=$2

ccefile=/home/jy344/Genomes/hg19/annotation/Ensemblv65/hg19_nh.chrom.sizes_tobins100_50_v65_merged_color.bed

cut -f 1-4 $ccefile | intersectBed -a $bedfile -b stdin -wo > ${bedfile/.bed/_intersect_cce.txt}
perl /home/jy344/Programs/Scripts/Age/calculate_age_by_cce_140411.pl ${bedfile/.bed/_intersect_cce.txt} $result
