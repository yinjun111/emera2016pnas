#Workflow for Figure 5, repeats and conservation

#For Figure 5A
#Merge MtAmn (Tetra, Gna, Vet) for rpts/conserv calculation
awk '{if ($4 == "Tetrapoda" || $4 == "Gnathostomata" || $4 == "Vertebrata") print $0 }' hg19_nh.chrom.sizes_tobins100_50_v65_merged_color.bed > hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmn.bed
awk '{if ($4 != "Tetrapoda" && $4 != "Gnathostomata" && $4 != "Vertebrata") print $0 }' hg19_nh.chrom.sizes_tobins100_50_v65_merged_color.bed > hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_NotMtAmn.bed
cut -f 1,2,3 hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmn.bed | mergeBed -i stdin | awk '{print $1"\t"$2"\t"$3"\tMtAmn"}' > hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmn_edi.bed
cut -f 1,2,3,4 hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_NotMtAmn.bed | sort - hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmn_edi.bed > hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmnMerged.bed


#Count cce coverage of whole genome
perl count_cce.pl #ccemerged_summary.txt
#Intersect cce with rpts
cut -f 1,2,3,4 /Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/MamCortex/genome/Human/cce_v65/hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmnMerged.bed | intersectBed -a stdin -b /Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/MamCortex/repbase/beds/hg19.rmsk.repeatmasker.bed -wo > cce2rpts_mtamnmerged.txt
perl count_cce_feature.pl cce2rpts_mtamnmerged.txt cce2rpts_mtamnmerged_summary.txt

#Count cce coverage of enhancers
perl count_enh_intersect_cce.pl #cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_lensummary.txt
perl count_enh_intersect_rpts.pl #cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_vs_rpts_lensummary.txt

#plot
Figure5A_rpts_summary.R



#For Figure 5C/D
#Intersect enhancers with cce, then repeats
for file in cortex_ac_v5_all.bed cortex_ac_v5_all_bypeak_enh_140723.bed;
do
	intersectBed -a $file -b hg19_nh.chrom.sizes_tobins100_50_v65_merged_color_MtAmnMerged.bed -wo > ${file/.bed/_intersect_ccemerged.txt}
done
#by peak shared enhancer
perl get_intersected_bed_141212.pl cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged.txt cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged.bed
#Intersect with repeats
intersectBed -a cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged.bed -b ~/Genomes/rmsk/hg19.rmsk.repeatmasker.bed -wo > cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_vs_rpts.txt
perl get_intersected_bed_141212.pl cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_vs_rpts.txt cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_vs_rpts.bed


#Intersect the enhancers/enhancers+cce/enhancers+cce+rpts with plancental phastcons element 
for file in cortex_ac_v5_all_*.bed
do
	intersectBed -a $file -b phastConsElements46wayPlacental.bed -wo > ${file/.bed/_intersect_phastConsElements46wayPlacental.txt}
done

#Summarize for phastcons ele #No intersection with PhastCons Ele will be marked as NA then treated as 0
perl summarize_phastconsel_for_cce_150223.pl cortex_ac_v5_all_bypeak_enh_140723.bed phastConsElements46wayPlacental

#plot
Figure5C_plot_conservation_El_for_CS23.R

