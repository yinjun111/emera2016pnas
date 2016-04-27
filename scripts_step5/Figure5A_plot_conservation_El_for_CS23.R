#

setwd("/Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/MamCortex/Repeats/Conservation_byEl")

data<-read.table("cortex_ac_v5_all_bypeak_enh_140723_CS23_phastconsel_phastConsElements46wayPlacental_summary_forccenrpts.txt",header=T,row.names=1,sep="\t")

clades<-c("Eutheria","Theria","Mammalia","Amniota","MtAmn")

#First plot by base

#all enh
data[is.na(data[,3]),3]=0
data.bybase.enh<-split(data[,3]/data[,2],data[,1])

#CCE #No PhastconsEl is treated as 0
data[is.na(data[,5]),5]=0
data.bybase.cce<-split(data[!is.na(data[,4]),5]/data[!is.na(data[,4]),4],data[!is.na(data[,4]),1])

#rpt
data[is.na(data[,7]),7]=0
data.bybase.rpt<-split(data[!is.na(data[,6]),7]/data[!is.na(data[,6]),6],data[!is.na(data[,6]),1])

#Figure 5C/D
#boxplot for percentage of no. of bases 
par(mfrow=c(1,3))
boxplot(data.bybase.enh[clades],ylab="Percentage of Bases included in the PhastCons Element",main="Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.bybase.cce[clades],ylab="Percentage of Bases included in the PhastCons Element",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.bybase.rpt[clades],ylab="Percentage of Bases included in the PhastCons Element",main="Repeats in the CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))

#Currently Used
par(cex=0.6)
boxplot(data.bybase.cce[clades],ylab="Percentage of Bases included in the PhastCons Element",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))



#Then by num
#---------------
data.bynum.enh<-rbind(table(data[,1])[clades],table(data[data[,3]!=0,1])[clades])
rownames(data.bynum.enh)<-c("Total No.","With PCE")

#CCE 
data.bynum.cce<-rbind(table(data[data[,4]!=0,1])[clades],table(data[data[,5]!=0,1])[clades])
rownames(data.bynum.cce)<-c("Total No.","With PCE")


#rpt
data.bynum.rpt<-rbind(table(data[data[,6]!=0,1])[clades],table(data[data[,7]!=0,1])[clades])
rownames(data.bynum.rpt)<-c("Total No.","With PCE")
 
#old
#---------------
par(mfrow=c(1,3))
coord<-barplot(rbind(data.bynum.enh[2,],data.bynum.enh[1,]-data.bynum.enh[2,]),col=c("blue","white"),axes=F,ylab="No. of features",main="Enhancers",ylim=c(0,3000))
axis(1,at=coord,labels=colnames(data.bynum.enh))
axis(2,at=seq(0,3000,500),labels=seq(0,3000,500))
axis(4,at=seq(0,1000,200),labels=seq(0,1000,200)/10,col="red",col.ticks="red")
lines(cbind(coord,data.bynum.enh[2,]/data.bynum.enh[1,]*1000),col="red")
box()
legend("topright",fill=c("white","blue"),legend=c("No. of features with PhastConsEl","No. of features without PhastConsEl"))

#cce
coord<-barplot(rbind(data.bynum.cce[2,],data.bynum.cce[1,]-data.bynum.cce[2,]),col=c("blue","white"),axes=F,ylab="No. of features",main="CladeConservedElement in the Enhancers",ylim=c(0,3000))
axis(1,at=coord,labels=colnames(data.bynum.cce))
axis(2,at=seq(0,3000,500),labels=seq(0,3000,500))
axis(4,at=seq(0,1000,200),labels=seq(0,1000,200)/10,col="red",col.ticks="red")
lines(cbind(coord,data.bynum.cce[2,]/data.bynum.cce[1,]*1000),col="red")
box()
legend("topright",fill=c("white","blue"),legend=c("No. of features with PhastConsEl","No. of features without PhastConsEl"))


#cce
coord<-barplot(rbind(data.bynum.rpt[2,],data.bynum.rpt[1,]-data.bynum.rpt[2,]),col=c("blue","white"),axes=F,ylab="No. of features",main="Repeats in the CladeConservedElement in the Enhancers",ylim=c(0,3000))
axis(1,at=coord,labels=colnames(data.bynum.rpt))
axis(2,at=seq(0,3000,500),labels=seq(0,3000,500))
axis(4,at=seq(0,1000,200),labels=seq(0,1000,200)/10,col="red",col.ticks="red")
lines(cbind(coord,data.bynum.rpt[2,]/data.bynum.rpt[1,]*1000),col="red")
box()
legend("topright",fill=c("white","blue"),legend=c("No. of features with PhastConsEl","No. of features without PhastConsEl"))


#new
#---------------

par(cex=0.6,mfrow=c(1,3))
coord<-barplot(data.bynum.enh[1,],col="white",axes=F,axisnames =F, ylab="No. of features",main="Enhancers",ylim=c(0,4000))
barplot(data.bynum.enh[2,],add=T,col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"),axes=F,axisnames =F)
axis(1,at=coord,labels=colnames(data.bynum.enh))
axis(2,at=seq(0,3000,500),labels=seq(0,3000,500))
axis(4,at=seq(0,3000,600)+200,labels=seq(0,100,20),col="red",col.ticks="red")
lines(cbind(coord,data.bynum.enh[2,]/data.bynum.enh[1,]*3000+200),col="red")
box()
legend("topright",fill=c("white","blue"),legend=c("No. of features without PhastConsEl","No. of features with PhastConsEl"))

par(cex=0.6)
coord<-barplot(data.bynum.cce[1,],col="white",axes=F,axisnames =F, ylab="No. of features",main="CladeConservedElement in the Enhancers",ylim=c(0,4000))
barplot(data.bynum.cce[2,],add=T,col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"),axes=F,axisnames =F)
axis(1,at=coord,labels=colnames(data.bynum.cce))
axis(2,at=seq(0,3000,500),labels=seq(0,3000,500))
axis(4,at=seq(0,3000,600)+200,labels=seq(0,100,20),col="red",col.ticks="red")
lines(cbind(coord,data.bynum.cce[2,]/data.bynum.cce[1,]*3000+200),col="red")
box()
legend("topright",fill=c("white","blue"),legend=c("No. of features without PhastConsEl","No. of features with PhastConsEl"))

par(cex=0.6)
coord<-barplot(data.bynum.rpt[1,],col="white",axes=F,axisnames =F, ylab="No. of features",main="Repeats in the CladeConservedElement in the Enhancers",ylim=c(0,4000))
barplot(data.bynum.rpt[2,],add=T,col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"),axes=F,axisnames =F)
axis(1,at=coord,labels=colnames(data.bynum.rpt))
axis(2,at=seq(0,3000,500),labels=seq(0,3000,500))
axis(4,at=seq(0,3000,600)+200,labels=seq(0,100,20),col="red",col.ticks="red")
lines(cbind(coord,data.bynum.rpt[2,]/data.bynum.rpt[1,]*3000+200),col="red")
box()
legend("topright",fill=c("white","blue"),legend=c("No. of features without PhastConsEl","No. of features with PhastConsEl"))




#Currently used!!
par(mfrow=c(1,2),cex=0.6)

#A
boxplot(data.bybase.cce[clades],ylab="Percentage of Bases included in the PhastCons Element",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))

#B
coord<-barplot(data.bynum.cce[1,],col="white",axes=F,axisnames =F, ylab="No. of features",main="CladeConservedElement in the Enhancers",ylim=c(0,4000))
barplot(data.bynum.cce[2,],add=T,col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"),axes=F,axisnames =F)
axis(1,at=coord,labels=colnames(data.bynum.cce))
axis(2,at=seq(0,3000,500),labels=seq(0,3000,500))
axis(4,at=seq(0,3000,600)+200,labels=seq(0,100,20),col="red",col.ticks="red")
lines(cbind(coord,data.bynum.cce[2,]/data.bynum.cce[1,]*3000+200),col="red")
box()
legend("topright",fill=c("#E8118A","white"),legend=c("No. of features with PhastConsEl","No. of features without PhastConsEl"))




