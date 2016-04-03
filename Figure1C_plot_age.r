#Figure 1C
#Pie chart

setwd("/Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/Repeat/ChIPseq/V5/Human/age")

data.ac<-read.table("../peaks_anno/cortex_ac_v5_all_annotation_combo_140723.txt",header=T,row.names=1,as.is=T,sep="\t")
anno.ac<-read.table("../DM/bed_annotation_ac.txt",header=F,row.names=NULL,as.is=T,sep="\t")

samples<-c("CS16","CS23","F2")
multiz.clades<-c("Human","Ape","Primate","Eutheria","Theria","Mammalia","Amniota","Tetrapoda","Gnathostomata","Vertebrata","NONE")
multiz.specs<-c("Human","Chimp","Orangutan","Rhesus","Marmoset","Mouse","Rat","Guinea Pig","Rabbit","Cow","Horse","Dog","Elephant","Opossum","Platypus","Chicken","Zebra finch","Lizard","Frog","Tetraodon","Fugu","Stickleback","Medaka","Zebrafish","Lamprey")


#Split the data by samples
data.ac.samples.clades<-list()
for(num in 1:nrow(anno.ac)) {
	data.ac.samples.clades[[anno.ac[num,3]]]<-data.ac[data.ac[,1]=="Enhancer" & grepl(anno.ac[num,2],data.ac[,5]) & grepl(anno.ac[num,3],data.ac[,4]),6]
}

#table
data.ac.samples.clades.table<-matrix(0,length(samples),length(multiz.clades))
for(sample.num in 1:length(samples)) {
	for( clade.num in 1:length(multiz.clades)) {
		data.ac.samples.clades.table[sample.num,clade.num]<-sum(data.ac.samples.clades[[samples[sample.num]]]==multiz.clades[clade.num])
	}
}
colnames(data.ac.samples.clades.table)<-multiz.clades
rownames(data.ac.samples.clades.table)<-samples

write.table(data.ac.samples.clades.table,file="cortex_ac_v5_all_clades_summary_140724.txt",sep="\t",quote=F,col.names=NA)



#Combine data
data.ac.samples.clades.table.combined<-matrix(0,length(samples),5)
data.ac.samples.clades.table.combined[,1:4]<-data.ac.samples.clades.table[,4:7]
data.ac.samples.clades.table.combined[,5]<-data.ac.samples.clades.table[,8]+data.ac.samples.clades.table[,9]+data.ac.samples.clades.table[,10]
rownames(data.ac.samples.clades.table.combined)<-samples
colnames(data.ac.samples.clades.table.combined)<-c(colnames(data.ac.samples.clades.table)[4:7],"MtAmiota")


#plot ages of the peaks
par(las=2,cex=0.8)
boxplot(data.samples.ages[samples],axes=F)
box()
axis(1,at=1:length(samples),labels=samples)
axis(2,at=1:25,labels=multiz.specs)

#plot clades
plot_pie<-function(anno) {
	par(mfrow=c(2,3))
	for(num in 1:nrow(anno)) {
		pie(as.numeric(anno[num,]),main=rownames(anno)[num],radius=1,labels =colnames(anno))
	}
}	
plot_pie(data.samples.clades.table.combined)
