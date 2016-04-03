#

setwd("/Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/MamCortex/Repeats/Conservation_byEl_byScore")

data<-read.table("cortex_ac_v5_all_bypeak_enh_140723_CS23_phastconsel_phastConsElements46wayPlacental_summary_forccenrpts_withlodnonorm_150205.txt",header=T,row.names=1,sep="\t")

clades<-c("Eutheria","Theria","Mammalia","Amniota","MtAmn")


#Figure 6B2
#5, total length
#---------------

#for enh without el, score is 0 or (no score, current) ?

#all enh
data[is.na(data[,3]),3]=0
data.byellen.enh<-split(data[,3],data[,1])

#CCE 
data[is.na(data[,7]),7]=0
data.byellen.cce<-split(data[,7],data[,1])

#rpt
data[is.na(data[,11]),11]=0
data.byellen.rpt<-split(data[,11],data[,1])


#old 
par(cex=0.6,mfrow=c(1,3))
boxplot(data.byellen.enh[clades],ylim=c(0,4000),ylab="Total Len of PhastconsEl",main="Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byellen.cce[clades],ylim=c(0,4000),ylab="Total Len of PhastconsEl",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byellen.rpt[clades],ylim=c(0,4000),ylab="Total Len of PhastconsEl",main="Repeats in the CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))

#current
par(cex=0.6,mfrow=c(1,2))
boxplot(list("Oldest part"=data.byellen.cce[["Mammalia"]],"Whole enhancer"=data.byellen.enh[["Mammalia"]]),ylim=c(0,800),ylab="Total bases in placental PhastCons Element",main="Mammalia",col=c("#F78E39","#F78E39"),outline=F)
boxplot(list("Oldest part"=data.byellen.cce[["MtAmn"]],"Whole enhancer"=data.byellen.enh[["MtAmn"]]),ylim=c(0,2100),ylab="Total bases in placental PhastCons Element",main="MtAmn",col=c("#C4E0AE","#C4E0AE"),outline=F)




#First plot by base
#---------------
#Fraction of bases coverd by El
#all enh
data[is.na(data[,3]),3]=0
data.bybase.enh<-split(data[,3]/data[,2],data[,1])

#CCE ##No PhastconsEl is treated as 0
data[is.na(data[,7]),7]=0
data.bybase.cce<-split(data[!is.na(data[,6]),7]/data[!is.na(data[,6]),6],data[!is.na(data[,6]),1])

#rpt
data[is.na(data[,11]),11]=0
data.bybase.rpt<-split(data[!is.na(data[,10]),11]/data[!is.na(data[,10]),10],data[!is.na(data[,10]),1])

#boxplot for percentage of no. of bases 
par(mfrow=c(1,3))
boxplot(data.bybase.enh[clades],ylab="Percentage of Bases included in the PhastCons Element",main="Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.bybase.cce[clades],ylab="Percentage of Bases included in the PhastCons Element",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.bybase.rpt[clades],ylab="Percentage of Bases included in the PhastCons Element",main="Repeats in the CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))

#Currently Used
par(cex=0.6)
boxplot(data.bybase.cce[clades],ylab="Percentage of Bases included in the PhastCons Element",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))

#wilcox test
data.bybase.cce.p<-c()
for (clade in clades) {
	data.bybase.cce.p<-c(data.bybase.cce.p,wilcox.test(data.bybase.cce[["Eutheria"]],data.bybase.cce[[clade]])$p.value)
}

#Then by num of enhancers with el
#---------------
data.bynum.enh<-rbind(table(data[,1])[clades],table(data[data[,3]!=0,1])[clades])
rownames(data.bynum.enh)<-c("Total No.","With PCE")

#CCE 
data.bynum.cce<-rbind(table(data[data[,6]!=0,1])[clades],table(data[data[,7]!=0,1])[clades])
rownames(data.bynum.cce)<-c("Total No.","With PCE")


#rpt
data.bynum.rpt<-rbind(table(data[data[,10]!=0,1])[clades],table(data[data[,11]!=0,1])[clades])
rownames(data.bynum.rpt)<-c("Total No.","With PCE")
 

#
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


#3. by highest score of el in enhancers
#---------------

#for enh without el, score is 0 or (no score, current) ?

#all enh
data.byhscore.enh<-split(data[!is.na(data[,5]),5],data[!is.na(data[,5]),1])
data.byhscore.enh2<-split(data[,5],data[,1])

#CCE 
data.byhscore.cce<-split(data[!is.na(data[,7]),9],data[!is.na(data[,7]),1])
data.byhscore.cce2<-split(data[,9],data[,1])

#rpt
data.byhscore.rpt<-split(data[!is.na(data[,11]),13],data[!is.na(data[,11]),1])
data.byhscore.rpt2<-split(data[,13],data[,1])

#
par(cex=0.6,mfrow=c(1,3))
boxplot(data.byhscore.enh[clades],ylim=c(0,300),ylab="Best length normalized PhastconsEl score",main="Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byhscore.cce[clades],ylim=c(0,300),ylab="Best length normalized PhastconsEl score",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byhscore.rpt[clades],ylim=c(0,300),ylab="Best length normalized PhastconsEl score",main="Repeats in the CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))



#3B. by highest score of el in enhancers
#---------------

#for enh without el, score is 0 or (no score, current) ?

#all enh
data.byhscore.enh<-split(data[!is.na(data[,5]),5],data[!is.na(data[,5]),1])
data.byhscore.enh2<-split(data[,5],data[,1])

#CCE 
data.byhscore.cce<-split(data[!is.na(data[,7]),9],data[!is.na(data[,7]),1])
data.byhscore.cce2<-split(data[,9],data[,1])

#rpt
data.byhscore.rpt<-split(data[!is.na(data[,11]),13],data[!is.na(data[,11]),1])
data.byhscore.rpt2<-split(data[,13],data[,1])

#
par(cex=0.6,mfrow=c(1,3))
boxplot(data.byhscore.enh[clades],ylab="Best length normalized PhastconsEl score",main="Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byhscore.cce[clades],ylab="Best length normalized PhastconsEl score",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byhscore.rpt[clades],ylab="Best length normalized PhastconsEl score",main="Repeats in the CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))



par(cex=0.6,mfrow=c(1,3))
boxplot(data.byhscore.enh[clades],ylab="Best PhastconsEl score",main="Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byhscore.cce[clades],ylab="Best PhastconsEl score",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byhscore.rpt[clades],ylab="Best PhastconsEl score",main="Repeats in the CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))

#wilcox test
data.byhscore.cce.p<-c()
for (clade in clades) {
	data.byhscore.cce.p<-c(data.byhscore.cce.p,wilcox.test(data.byhscore.cce[["Eutheria"]],data.byhscore.cce[[clade]])$p.value)
}


#4. by number of els in enhancers
#---------------


#for enh without el, score is 0 or (no score, current) ?

#all enh
data[is.na(data[,4]),4]=0
data.byelnum.enh<-split(data[,4],data[,1])

#CCE 
data[is.na(data[,8]),8]=0
data.byelnum.cce<-split(data[,8],data[,1])

#rpt
data[is.na(data[,12]),12]=0
data.byelnum.rpt<-split(data[,12],data[,1])


#
par(cex=0.6,mfrow=c(1,3))
boxplot(data.byelnum.enh[clades],ylim=c(0,60),ylab="No. of PhastconsEl",main="Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byelnum.cce[clades],ylim=c(0,60),ylab="No. of PhastconsEl",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byelnum.rpt[clades],ylim=c(0,60),ylab="No. of PhastconsEl",main="Repeats in the CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))



#5, total length
#---------------

#for enh without el, score is 0 or (no score, current) ?

#all enh
data[is.na(data[,3]),3]=0
data.byellen.enh<-split(data[,3],data[,1])

#CCE 
data[is.na(data[,7]),7]=0
data.byellen.cce<-split(data[,7],data[,1])

#rpt
data[is.na(data[,11]),11]=0
data.byellen.rpt<-split(data[,11],data[,1])



par(cex=0.6,mfrow=c(1,3))
boxplot(data.byellen.enh[clades],ylim=c(0,4000),ylab="Total Len of PhastconsEl",main="Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byellen.cce[clades],ylim=c(0,4000),ylab="Total Len of PhastconsEl",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byellen.rpt[clades],ylim=c(0,4000),ylab="Total Len of PhastconsEl",main="Repeats in the CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))

#Avg len
data[data[,4]==0,4]=0.1
data.byelavglen.enh<-split(data[,3]/data[,4],data[,1])

#CCE 
data[data[,8]==0,8]=0.1

data.byelavglen.cce<-split(data[,7]/data[,8],data[,1])

#rpt
data[data[,12]==0,12]=0.1
data.byelavglen.rpt<-split(data[,11]/data[,12],data[,1])


#
par(cex=0.6,mfrow=c(1,3))
boxplot(data.byelavglen.enh[clades],ylim=c(0,500),ylab="Avg Len of PhastconsEl",main="Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byelavglen.cce[clades],ylim=c(0,500),ylab="Avg Len of PhastconsEl",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
boxplot(data.byelavglen.rpt[clades],ylim=c(0,500),ylab="Avg Len of PhastconsEl",main="Repeats in the CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))




#-------------------



#Currently used!!
par(mfrow=c(1,2),cex=0.8)

#A
boxplot(data.bybase.cce[clades],ylab="Percentage of Bases included in the PhastCons Element",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
#points((1:5)[data.bybase.cce.p<0.01],rep(1,sum(data.bybase.cce.p<0.01)),pch="*",col="red")


#B
boxplot(data.byhscore.cce[clades],ylim=c(0,2700),ylab="Best PhastconsEl LOD score",main="CladeConservedElement in the Enhancers",col=c("#E8118A","#2E3192","#F58D39","#00AEEF","#009444"))
#points((1:5)[data.byhscore.cce.p<0.01],rep(2700,sum(data.byhscore.cce.p<0.01)),pch="*",col="red")





