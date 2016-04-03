#Workflow for Figure 6, conservation


#For Figure 6A/B

#Get aggregated bed file
cat hg19_nh.chrom.sizes_tobins100_50_v65_Eutheria_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_Theria_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_Mammalia_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_Amniota_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmn_edi.bed | sort -k1,1 -k2,2n | mergeBed -i stdin | awk '{print $1"\t"$2"\t"$3"\tEutheriaNOlder"}' > hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_EutheriaNOlder_merged.bed
cat hg19_nh.chrom.sizes_tobins100_50_v65_Theria_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_Mammalia_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_Amniota_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmn_edi.bed | sort -k1,1 -k2,2n | mergeBed -i stdin | awk '{print $1"\t"$2"\t"$3"\tTheriaNOlder"}' > hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_TheriaNOlder_merged.bed
cat hg19_nh.chrom.sizes_tobins100_50_v65_Mammalia_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_Amniota_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmn_edi.bed | sort -k1,1 -k2,2n | mergeBed -i stdin | awk '{print $1"\t"$2"\t"$3"\tMammaliaNOlder"}' > hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MammaliaNOlder_merged.bed
cat hg19_nh.chrom.sizes_tobins100_50_v65_Amniota_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmn_edi.bed | sort -k1,1 -k2,2n | mergeBed -i stdin | awk '{print $1"\t"$2"\t"$3"\tAmiotaNOlder"}' > hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_AmiotaNOlder_merged.bed
cat hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_EutheriaNOlder_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_TheriaNOlder_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MammaliaNOlder_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_AmiotaNOlder_merged.bed hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmn_edi.bed > hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_AggreMerged.bed

#intersect age segments with cce
intersectBed -a cortex_ac_v5_all_bypeak_enh_140723.bed -b /Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/MamCortex/genome/Human/cce_v65/hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_AggreMerged.bed -wo > cortex_ac_v5_all_bypeak_enh_140723_intersect_cceaggremerged.txt
perl get_intersected_bed_141212.pl cortex_ac_v5_all_bypeak_enh_140723_intersect_cceaggremerged.txt cortex_ac_v5_all_bypeak_enh_140723_intersect_cceaggremerged.bed

#Get the annotation suitable for the following analysis
perl get_bypeak_anno.pl #cortex_ac_v5_all_annotation_clade_all_forcount_140723_selbypeak.txt

#Summarize the length of different enhancers
perl summarize_length_by_cce_150609.pl


#Plot Figure 6A/B
Figure6AB_plot_len_byCCE.R



#Figure 6B2, total bases
perl summarize_phastconsel_for_cce_150225_lod.pl cortex_ac_v5_all_bypeak_enh_140723.bed phastConsElements46wayPlacental
#cortex_ac_v5_all_bypeak_enh_140723_CS23_phastconsel_phastConsElements46wayPlacental_summary_forccenrpts_withlodnonorm_150205.txt

#Plot Figure 6B2
Figure6B2_plot_conservation_El_for_CS23.R
#cortex_ac_v5_all_bypeak_enh_140723_CS23_phastconsel_phastConsElements46wayPlacental_summary_forccenrpts_150203_byTotalLen.pdf

