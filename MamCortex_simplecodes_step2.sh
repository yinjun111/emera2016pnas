#Workflow to age the whole genome, defined shuffle background, and get tracks for UCSC genome browser

#Get the non-exonic regions in the genome
subtractBed -a hg19_nh.chrom.sizes_tobins100_50.bed -b ~/Genomes/hg19/annotation/Ensemblv65/Homo_sapiens.gencode.v10.compositeexon.bed > hg19_nh.chrom.sizes_tobins100_50_v65_forSubmit.bed

#Split the bed files, run the process in multiple threads
split -l 1000000 -d hg19_nh.chrom.sizes_tobins100_50_forSubmit.bed hg19_nh.chrom.sizes_tobins100_50_forSubmit.bed_p
#make pbs
for file in *bed_p*;do sh /home/jy344/Programs/Scripts/Age/calculate_age_131121.sh $file;done;
cat hg19_nh.chrom.sizes_tobins100_50_forSubmit_submit.pbs_p0* hg19_nh.chrom.sizes_tobins100_50_forSubmit_submit.pbs_p1* > submit1.pbs
cat hg19_nh.chrom.sizes_tobins100_50_forSubmit_submit.pbs_p2* hg19_nh.chrom.sizes_tobins100_50_forSubmit_submit.pbs_p3* > submit2.pbs
cat hg19_nh.chrom.sizes_tobins100_50_forSubmit_submit.pbs_p4* hg19_nh.chrom.sizes_tobins100_50_forSubmit_submit.pbs_p5* hg19_nh.chrom.sizes_tobins100_50_forSubmit_submit.pbs_p6* > submit3.pbs


#process the results
split -l 1130000 -d hg19_nh.chrom.sizes_tobins100_50_ageresults.txt hg19_nh.chrom.sizes_tobins100_50_ageresults.txt_p
for file in hg19_nh.chrom.sizes_tobins100_50_ageresults.txt_p*;do echo "cd /home/jy344/Scratch/Projects/Repeats/genome/age10k;perl /home/jy344/Programs/Scripts/Age/process_ageresults_131125.pl $file ${file/.txt/_processed.txt} 0.5" >> process.pbs;done;
cat hg19_nh.chrom.sizes_tobins100_50_ageresults_processed.txt_p* > hg19_nh.chrom.sizes_tobins100_50_ageresults_all_processed.txt

#summarize the results
split -l 570000 -d hg19_nh.chrom.sizes_tobins100_50_ageresults_all_processed.txt hg19_nh.chrom.sizes_tobins100_50_ageresults_all_processed.txt_p
if [ -r summary.pbs ];then rm summary.pbs;fi;
for file in hg19_nh.chrom.sizes_tobins100_50_ageresults_all_processed.txt_p*;do echo "cd /home/jy344/Scratch/Projects/Repeats/genome/age10k;perl /home/jy344/Programs/Scripts/Age/summarize_ageresults_131212.pl $file ${file/_all_processed.txt/_summary.txt}" >> summary.pbs;done;
cat hg19_nh.chrom.sizes_tobins100_50_ageresults_summary.txt_* > hg19_nh.chrom.sizes_tobins100_50_ageresults_all_summary.txt

#find lost windows, e.g. windows without multiz alignment
perl ~/Programs/Scripts/Age/findlost_ageresults_131125.pl hg19_nh.chrom.sizes_tobins100_50_ageresults_all_summary.txt hg19_nh.chrom.sizes_tobins100_50_forSubmit.bed hg19_nh.chrom.sizes_tobins100_50_ageresults_all_summary_edi.t_all_summary_edi.txt

#Merge age analysis into elements
mkdir merged notmerged
perl /home/jy344/Programs/Scripts/Age/split_ageresults_byclade_131213.pl ../age10k_v65/hg19_nh.chrom.sizes_tobins100_50_ageresults_v65_all_all_summary_edi.txt ../age10k_v65/hg19_nh.chrom.sizes_tobins100_50_v65_forSubmit.bed notmerged/hg19_nh.chrom.sizes_tobins100_50_v65
for file in notmerged/*.bed;
    do filename=`echo $file | cut -d"/" -f2`;cut -f 1,2,3 $file | mergeBed -i stdin | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > merged/${filename/.bed/_merged.bed};
done;

#Summarize the elements, and color the beds
perl convert_cce_ucsc_colorbed.pl #hg19_nh.chrom.sizes_tobins100_50_v65_merged_color.bed




#Convert age results into UCSC tracks
#bed 2 bedgraph
if [ -r bg.pbs ];then rm bg.pbs;fi;
for file in notmerged/*.bed;do filename=`echo $file | cut -d"/" -f2`;echo "cd /home/jy344/Scratch/Projects/Repeats/genome/cce_v65;perl /home/jy344/Programs/Scripts/Age/convert_cce_ucsc_colorbedgraph_131209.pl $file bedgraph/${filename/.bed/_color.bg}" >> bg.pbs;done;

cat bedgraph/hg19_nh.chrom.sizes_tobins100_50_* | sort -k1,1 -k2,2n > hg19_nh.chrom.sizes_tobins100_50_all_color_sorted.bedgraph

bedGraphToBigWig hg19_nh.chrom.sizes_tobins100_50_all_color_sorted.bedgraph /home/jy344/Genomes/hg19/annotation/hg19_nh.chrom.sizes hg19_nh.chrom.sizes_tobins100_50_all_color_sorted.bw

#merge beds into color bed
perl convert_cce_ucsc_colorbed.pl
#bed to bigbed for ucsc
bedToBigBed hg19_nh.chrom.sizes_tobins100_50_v65_merged_color.bed ../v67/hg19_nh.chrom.sizes hg19_nh.chrom.sizes_tobins100_50_v65_merged_color.bb





