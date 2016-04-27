#
setwd("/Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/MamCortex/Conservation/Bases")

#genome background
genomeccebk<-read.table("ccemerged_summary.txt",header=T,row.names=1,sep="\t")
genomerptsbk<-read.table("cce2rpts_mtamnmerged_summary.txt",header=T,row.names=1,sep="\t")

genomeccebk.edi<-genomeccebk[c(4,2,11,3,12,5,1,6),]
genomerptsbk.edi<-genomerptsbk[c(4,2,11,3,12,5,1,6),]

#shared enhancer
sharedenhcce<-read.table("cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_lensummary.txt",header=T,row.names=1,sep="\t")
sharedenhrpts<-read.table("cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_vs_rpts_lensummary.txt",header=T,row.names=1,sep="\t")

#plot rpts contribution
#genoem background
genome.bar<-cbind(genomerptsbk.edi[,2],genomeccebk.edi[,2]-genomerptsbk.edi[,2])[2:8,]
rownames(genome.bar)<-colnames(sharedenhcce)[2:8]
colnames(genome.bar)<-c("RPTS","Not RPTS")

#CS23 shared enh
sharedenh.cs23.bar<-cbind(unlist(sharedenhrpts[2,]),unlist(sharedenhcce[2,]-sharedenhrpts[2,]))[2:8,]
colnames(sharedenh.cs23.bar)<-c("RPTS","Not RPTS")

#New!!
#------------
par(mfrow=c(1,2),cex=0.6,las=3)

coord<-barplot(genome.bar[,1]+genome.bar[,2],col="white",axes=F,axisnames=F,ylab="Billion Bases",ylim=c(0,1700000000),main="Whole Genome")
barplot(genome.bar[,1],col=c("grey","grey","#E8118A","#2E3192","#F58D39","#00AEEF","#009444"),axes=F,axisnames=F,add=T)

axis(1,at=coord,labels=rownames(genome.bar))
axis(2,at=seq(0,1500000000,500000000),labels=seq(0,1.5,0.5))
#axis(4,at=seq(0,1000000000,200000000),labels=seq(0,100,20),col="blue",col.ticks="blue")
#lines(cbind(coord,genome.bar[,1]/(genome.bar[,1]+genome.bar[,2])*1000000000),col="blue")
box()
legend("topright",fill=c("#E8118A","white"),legend=c("Bases covered by rpts","Bases NOT covered by rpts"))

coord<-barplot(sharedenh.cs23.bar[,1]+sharedenh.cs23.bar[,2],col="white",axes=F,axisnames=F,ylab="Million Bases",ylim=c(0,15000000),main="CS23 Human&Mouse shared enhancers")
barplot(sharedenh.cs23.bar[,1],col=c("grey","grey","#E8118A","#2E3192","#F58D39","#00AEEF","#009444"),axes=F,axisnames=F,add=T)

axis(1,at=coord,labels=rownames(sharedenh.cs23.bar))
axis(2,at=seq(0,15000000,5000000),labels=seq(0,1.5,0.5))
#axis(4,at=seq(0,10000000,2000000),labels=seq(0,100,20),col="blue",col.ticks="blue")
#lines(cbind(coord,sharedenh.cs23.bar[,1]/(sharedenh.cs23.bar[,1]+sharedenh.cs23.bar[,2])*10000000),col="blue")
box()
legend("topright",fill=c("#E8118A","white"),legend=c("Bases covered by rpts","Bases NOT covered by rpts"))

