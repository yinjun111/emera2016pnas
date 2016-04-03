#From fasta, alignment, peak calling to peak annotation

#Merge bowtie alignments for Human F2F, F2O, Mouse e17.5F e17.5O
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment;rm e17.5_ac_rep1.bowtie;zcat /home/skr2/Scratch/Steve/bowties/SR010_BC12.bowtie.gz >> e17.5_ac_rep1.bowtie;cat /home/skr2/Scratch/Steve/SR022_BC7.bowtie >> e17.5_ac_rep1.bowtie
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment;zcat /home/skr2/Scratch/Steve/bowties/SR035_BC6.bowtie.gz /home/skr2/Scratch/Steve/bowties/SR024_BC6.bowtie.gz > e17.5_ac_rep2.bowtie
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment;zcat /home/skr2/Scratch/Steve/bowties/SR009-BC2.bowtie.gz /home/skr2/Scratch/Steve/bowties/SR027_BC5.bowtie.gz > e17.5_input_rep1.bowtie
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment;zcat /home/skr2/Scratch/Steve/bowties/SR034_BC2.bowtie.gz /home/skr2/Scratch/Steve/bowties/SR023_BC7.bowtie.gz > e17.5_input_rep2.bowtie
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment;zcat /home/skr2/Scratch/Steve/bowties/SR028_BC5.bowtie.gz /home/skr2/Scratch/Steve/bowties/SR033_BC4.bowtie.gz > F2_ac_rep1.bowtie
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment;zcat /home/skr2/Scratch/Steve/bowties/SR011-BC2.bowtie.gz /home/skr2/Scratch/Steve/bowties/SR033_BC12.bowtie.gz > F2_ac_rep2.bowtie
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment;rm F2_input_rep1.bowtie;zcat /home/skr2/Scratch/Steve/bowties/SR020_BC4.bowtie.gz >> F2_input_rep1.bowtie;cat /home/skr2/Scratch/Steve/SR023_BC2.bowtie >> F2_input_rep1.bowtie
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment;rm F2_input_rep2.bowtie;zcat /home/skr2/Scratch/Steve/bowties/SR011-BC12.bowtie.gz >> F2_input_rep2.bowtie;cat /home/skr2/Scratch/Steve/SR012.1_BC4.bowtie >> F2_input_rep2.bowtie

#call peaks only for Human F2 and Mouse e17.5. The other peaks are from Steve
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/peaks; perl /home/jy344/Programs/Scripts/Steve_Windower_bowtie_JLedit_JY_130701.pl /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment/F2_ac_rep1 /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment/F2_input_rep1 /home/jy344/Genomes/hg19/annotation/hg19_nh.chrom.sizes 0.00001 8
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/peaks; perl /home/jy344/Programs/Scripts/Steve_Windower_bowtie_JLedit_JY_130701.pl /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment/F2_ac_rep2 /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment/F2_input_rep2 /home/jy344/Genomes/hg19/annotation/hg19_nh.chrom.sizes 0.00001 8
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/peaks; perl /home/jy344/Programs/Scripts/Steve_Windower_bowtie_JLedit_JY_130701.pl /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment/e17.5_ac_rep1 /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment/e17.5_input_rep1 /home/jy344/Genomes/mm9/annotation/v67/mm9_nh.chrom.sizes 0.00001 8
cd /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/peaks; perl /home/jy344/Programs/Scripts/Steve_Windower_bowtie_JLedit_JY_130701.pl /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment/e17.5_ac_rep2 /home/jy344/Scratch/Projects/Repeats/ChIPseq/V5/alignment/e17.5_input_rep2 /home/jy344/Genomes/mm9/annotation/v67/mm9_nh.chrom.sizes 0.00001 8

cat F2*ac*bed |mergeBed -i stdin -n |awk '{if($4>1) print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >merge_F2_ac_overlap.bed
cat F2*me2*bed |mergeBed -i stdin -n |awk '{if($4>1) print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >merge_F2_me2_overlap.bed
cat e17.5*ac*bed |mergeBed -i stdin -n |awk '{if($4>1) print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >merge_e17_ac_overlap.bed
cat e17.5*me2*bed |mergeBed -i stdin -n |awk '{if($4>1) print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >merge_e17_me2_overlap.bed

cat merge_CS*overlap.bed merge_F*overlap.bed | grep -v "chrM" | grep -v "_" | grep -v "\." |  awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' | sort | uniq > cortex_v5_all.bed
cat merge_e11_overlap_mm.bed merge_e14_overlap_mm.bed merge_e17_ac_overlap.bed | grep -v "chrM" | grep -v "_" | grep -v "\." |  awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' | sort | uniq > cortex_ac_mm_v5_all.bed


#Intersect peaks 
#Reciprocal match the peaks from human to mouse & from mouse to human
sh ~/Programs/Scripts/Bed/reciprocal_intersect_130714.sh cortex_ac_v5_all.bed cortex_ac_mm_v5_all.bed /home/jy344/Genomes/alignments/hg19ToMm9.over.chain.gz /home/jy344/Genomes/alignments/mm9ToHg19.over.chain.gz cortex_v5_hs_mm_intersect_recip_0.txt 0

#Annotation of peaks by exon/promoter/enhancer, distance to genes, then by tissues
sh annobed_regions_v65_140423.sh ../peaks/cortex_ac_v5_all.bed cortex_ac_v5_all_annoregion.txt
sh annobed_v65_140423.sh ../peaks/cortex_ac_v5_all.bed cortex_ac_v5_all_annoez.txt
sh ~/Scripts/Bed/annobed_regions_general_140425.sh ../peaks/cortex_ac_mm_v5_all.bed ~/Projects/Repeat/genome/Mouse/v65/Mus_musculus.NCBIM37.65_edi_annotated.bed cortex_ac_mm_v5_all_annoregion.txt
sh ~/Scripts/Bed/annobed_ez_general_140423.sh  ../peaks/cortex_ac_mm_v5_all.bed ~/Projects/Repeat/genome/Mouse/v65/Mus_musculus.NCBIM37.65_edi_genemerged_edi.bed ~/Projects/Repeat/genome/Mouse/v65/Mus_musculus.NCBIM37.65_edi_promoter_1000_edi.bed cortex_ac_mm_v5_all_annoez.txt
perl dm_byPeak_cortex_ac.pl #get cortex_ac_v5_all_peakanno140722.txt
perl dm_byPeak_cortex_ac_mm.pl #get cortex_ac_mm_v5_all_peakanno_140722.txt

#Get Functional Conserved peaks by peaks overlap in human & mouse
perl dm_by_cross_spec.pl cortex_ac_v5_all_peakanno140722.txt cortex_ac_mm_v5_all_peakanno_140722.txt ../peaks_intersect/cortex_v5_hs_mm_intersect_recip_0.txt cortex_ac_v5_all_human_tomouse_p0_bothdm_byPeaks_140722.txt 3 inter

#Age analysis #calculate age by intersecting with CCE
sh ~/Programs/Scripts/Age/calculate_age_by_cce_140411.sh cortex_ac_v5_all.bed cortex_ac_v5_all_age.txt

#Lastly summarize the annotation
perl combo-anno_140723_ac.pl #get cortex_ac_v5_all_annotation_combo_140723.txt




