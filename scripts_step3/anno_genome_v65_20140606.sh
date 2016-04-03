#annotate the genome based on the GTF file

exonfile="gencode.v10.wholeGene.exonTranscript_edi.gtf"
gapfile="hg19_Gap.bed"
chrsize="hg19_nh.chrom.sizes"
promoterupstream=1000


##GTF 2 bed
#perl ~/Scripts/Format/gtf2bed_120626.pl $exonfile ${exonfile/.gtf/.bed}

#echo "Identify promoters!\n"
#perl ~/Scripts/Format/gtf2ucsc_120501.pl $exonfile ${exonfile/.gtf/_ucsc.txt}
#perl ~/Scripts/get_promoter_common_130404.pl --input ${exonfile/.gtf/_ucsc.txt} --upstream $promoterupstream --output ${exonfile/.gtf/_promoter_$promoterupstream.bed} --feature tx

#echo "Merge exons!\n"
#mkdir merged genes;for gene in `cut -f 4 ${exonfile/.gtf/.bed}  | sort | uniq`;do grep $gene ${exonfile/.gtf/.bed} > genes/$gene.bed;mergeBed -nms -i genes/$gene.bed > merged/$gene.bed;done;
#for file in merged/*.bed;do cat $file >> ${exonfile/.gtf/_genemerged.bed};done;

##final merged exon set #changed 140423
#perl ~/Scripts/Format/dedup_merged_exon.pl ${exonfile/.gtf/_genemerged.bed} ${exonfile/.gtf/_genemerged_dedup.bed}
##this caused some genes missing
##grep -v "_" ${exonfile/.gtf/_genemerged_dedup.bed} | grep -v "\." | awk '{if($1=="chrMT") print "chrM\t"$2"\t"$3"\t"$4; else print $0}' | sort -k 1,1 -k 2,2n > ${exonfile/.gtf/_genemerged_edi.bed}
cat ${exonfile/.gtf/_genemerged_dedup.bed} | awk '{if($1=="chrMT") print "chrM\t"$2"\t"$3"\t"$4; else print $0}' | sort -k 1,1 -k 2,2n > ${exonfile/.gtf/_genemerged_edi.bed}

#this step is not necessary sometimes..
##grep -v "_" ${exonfile/.gtf/_promoter_$promoterupstream.bed} | grep -v "\." > ${exonfile/.gtf/_promoter_$promoterupstream\_edi.bed}
cat ${exonfile/.gtf/_promoter_$promoterupstream.bed} > ${exonfile/.gtf/_promoter_$promoterupstream\_edi.bed}

#merged exon+promoter
cat ${exonfile/.gtf/_genemerged_edi.bed} ${exonfile/.gtf/_promoter_$promoterupstream\_edi.bed} | cut -f 1,2,3 | mergeBed -i stdin > ${exonfile/.gtf/_exonpromoter.bed}

echo "Get regions!\n"
#intron
perl ~/Scripts/Bed/get_gene_regions_140425.pl ${exonfile/.gtf/_genemerged_edi.bed} $chrsize ${exonfile/.gtf/_genemerged_edi_intron.bed} intron
subtractBed -a ${exonfile/.gtf/_genemerged_edi_intron.bed} -b ${exonfile/.gtf/_exonpromoter.bed} | subtractBed -a stdin -b $gapfile | mergeBed -i stdin > ${exonfile/.gtf/_genemerged_Intron.bed}

#10k
perl ~/Scripts/Bed/get_gene_regions_140425.pl ${exonfile/.gtf/_genemerged_edi.bed} $chrsize ${exonfile/.gtf/_genemerged_edi_10000.bed} 10000
#minus exon promoter intron
cat ${exonfile/.gtf/_genemerged_Intron.bed} ${exonfile/.gtf/_exonpromoter.bed} | awk '{print $1"\t"$2"\t"$3}'> ${exonfile/.gtf/_genemerged_no10k.bed}
subtractBed -a ${exonfile/.gtf/_genemerged_edi_10000.bed} -b ${exonfile/.gtf/_genemerged_no10k.bed} | subtractBed -a stdin -b $gapfile | mergeBed -i stdin > ${exonfile/.gtf/_genemerged_10k.bed}

#100k
perl ~/Scripts/Bed/get_gene_regions_140425.pl ${exonfile/.gtf/_genemerged_edi.bed} $chrsize ${exonfile/.gtf/_genemerged_edi_100000.bed} 100000
cat ${exonfile/.gtf/_genemerged_Intron.bed} ${exonfile/.gtf/_exonpromoter.bed} ${exonfile/.gtf/_genemerged_10k.bed}| awk '{print $1"\t"$2"\t"$3}'> ${exonfile/.gtf/_no100k.bed}
subtractBed -a ${exonfile/.gtf/_genemerged_edi_100000.bed} -b ${exonfile/.gtf/_no100k.bed} | subtractBed -a stdin -b $gapfile | mergeBed -i stdin > ${exonfile/.gtf/_genemerged_100k.bed}

#1-100k
subtractBed -a ${exonfile/.gtf/_genemerged_edi_100000.bed} -b ${exonfile/.gtf/_genemerged_no10k.bed} | subtractBed -a stdin -b $gapfile | mergeBed -i stdin > ${exonfile/.gtf/_genemerged_100ks.bed}

#mt100k
cat ${exonfile/.gtf/_no100k.bed} ${exonfile/.gtf/_genemerged_100k.bed}| awk '{print $1"\t"$2"\t"$3}' > ${exonfile/.gtf/_nomt100k.bed}
awk '{print $1"\t0\t"$2}' $chrsize > $chrsize.bed
subtractBed -a $chrsize.bed -b ${exonfile/.gtf/_nomt100k.bed} | subtractBed -a stdin -b $gapfile | mergeBed -i stdin > ${exonfile/.gtf/_genemerged_mt100k.bed}

echo "Final Step!\n"
#annotation
if [ -r ${exonfile/.gtf/_annotated.bed} ];
then
	rm ${exonfile/.gtf/_annotated.bed}
fi
awk '{print $1"\t"$2"\t"$3"\tExon"}' ${exonfile/.gtf/_genemerged_edi.bed} >> ${exonfile/.gtf/_annotated.bed}
awk '{print $1"\t"$2"\t"$3"\tPromoter"}' ${exonfile/.gtf/_promoter_$promoterupstream.bed} >> ${exonfile/.gtf/_annotated.bed}
awk '{print $1"\t"$2"\t"$3"\tIntron"}' ${exonfile/.gtf/_genemerged_Intron.bed} >> ${exonfile/.gtf/_annotated.bed}
awk '{print $1"\t"$2"\t"$3"\t10k"}' ${exonfile/.gtf/_genemerged_10k.bed} >> ${exonfile/.gtf/_annotated.bed}
awk '{print $1"\t"$2"\t"$3"\t100k"}'  ${exonfile/.gtf/_genemerged_100k.bed} >> ${exonfile/.gtf/_annotated.bed}
awk '{print $1"\t"$2"\t"$3"\tmt100k"}' ${exonfile/.gtf/_genemerged_mt100k.bed} >> ${exonfile/.gtf/_annotated.bed}
