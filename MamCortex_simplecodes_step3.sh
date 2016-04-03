#Workflow for permutation test, including shuffling, intersection and counting

#First prepare the background for shuffling
#Get the exon, promoter, intron, 10k, 100k, mt 100k region annotation
sh anno_genome_v65_20140606.sh
#Merge the distance based regions
for file in ../../v65/regions/gencode.v10.wholeGene.exonTranscript_edi_*.bed; do filename=`echo $file | awk 'BEGIN {FS="/"}; {print $NF}'`; mergeBed -i $file > regions_merged/$filename;done;
#Intersect the distance based regions with Age segments
perl prepare_cce_for_shuffling.pl
sh gen_region_forcce.sh
#Process the intersection file, give names to each distance-age segments
for file in intersect/*_intersect.txt;do filename=`echo $file | awk 'BEGIN {FS="/"}; {print $NF}'`;perl process_cce_intersect.pl $file  background/${filename/_intersect.txt/_shuffle_background.txt};done;
cat background/*_shuffle_background.txt > gencode_v10_cce_all_shuffle_background.txt
#Only need the distance-age segments for >eutheria and enhancer
grep -v "Exon" gencode_v10_cce_all_shuffle_background.txt | grep -v "Promoter" | grep -v "Human" | grep -v "Primate" | grep -v "Ape"  > gencode_v10_cce_all_shuffle_background_short.txt
#Give index for the distance-age segments for shuffling
perl preprocess_cce_for_shuffle.pl   #gencode_v10_cce_all_shuffle_background_forShuffle_140606.txt


#Get the bed files available (only use Eutherian and older) for shuffling
perl get_bed_for_shuffle_ac.pl
#Preprocesss the bed file, index by clade, region and length
#coord in the error.log may not get unique location when shuffled together
perl preprocess_bed_for_shuffle.pl cortex_ac_v5_all_bypeak_enh_140723.bed /Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/Repeat/ChIPseq/V5/Human/peaks_anno/cortex_ac_v5_all_annotation_combo_140723.txt /Users/junyin/Projects/Repeat/genome/Human/cce_v65/background/gencode_v10_cce_all_shuffle_background_forShuffle_140606.txt cortex_ac_v5_all_bypeak_enh_140723_anno_forshuffle.txt >& error.log


#ready to shuffle
perl gen_shuffle_workflow_cce_v1.pl
