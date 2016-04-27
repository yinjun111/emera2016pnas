#
setwd("/Volumes/Illumina1/Users1/jyin/WorkYale_2012/Projects/MamCortex/Conservation/Conservation_byCCE")


#By Placental PhastCons Element

data<-read.table("cortex_ac_v5_all_bypeak_enh_140723_intersect_ccemerged_length_bycceaggresummary.txt",header=T,row.names=1,sep="\t")
anno<-read.table("cortex_ac_v5_all_annotation_clade_all_forcount_140723_selbypeak.txt",header=T,row.names=1,sep="\t")

for(sample in c("CS16_embryonic_cortex","CS23_embryonic_cortex","F2_embryonic_cortex")) {

	pdf(paste("plot_aggregated_length_byCCE_",sample,".pdf",sep=""),width=12,height=3)
	par(las=2,cex=0.8,mar=c(8, 8, 4, 2),mfrow=c(1,5))
	
	for(num1 in 1:length(clades)) {
		clade1=clades[num1] #current clade
		if(clade1=="MtAmn") {
			data.sel<-data[grepl(sample,anno[,4]) & (grepl("Tetrapoda",anno[,5]) | grepl("Gnathostomata",anno[,5]) | grepl("Vertebrata",anno[,5]) ),]
		}
		else {
			data.sel<-data[grepl(sample,anno[,4]) & grepl(clade1,anno[,5]),]
		}
		
		plot(0,xlim=c(-2500,2500),ylim=c(0,7),pch="",main=clade1,axes=F,xlab="Aggregated Length",ylab="")
		
		for(num2 in 1:(num1+1)) {
			len.median<-median(data.sel[data.sel[,num2]>0,num2])
			len.quartiles<-quantile(data.sel[data.sel[,num2]>0,num2])
			rect(-len.median/2,num2-0.5,len.median/2,num2+0.5)
			#left
			#75 perc
			draw_horizontal_errbar(-len.median/2,num2,-(len.median/2+(len.quartiles[4]-len.median)/2),num2,0.2)
			#25 perc
			draw_horizontal_errbar(-len.median/2,num2,-(len.median/2-(len.median-len.quartiles[2])/2),num2,0.2)			
			#right
			#75 perc
			draw_horizontal_errbar(len.median/2,num2,len.median/2+(len.quartiles[4]-len.median)/2,num2,0.2)
			#25 perc
			draw_horizontal_errbar(len.median/2,num2,len.median/2-(len.median-len.quartiles[2])/2,num2,0.2)				
		}
		
		axis(1,at=seq(-2500,2500,500),labels=seq(-2500,2500,500))
		axis(2,at=1:6,labels=c("Whole Enhancer",clades))
		box()
		
	}
	dev.off()
	
}


#only plot mammalia and amiota
draw_horizontal_errbar<-function(x1,y1,x2,y2,cap) {
	segments(x1,y1,x2,y2)
	segments(x2,y2-cap/2,x2,y2+cap/2)
}

draw_vertical_errbar<-function(x1,y1,x2,y2,cap) {
	segments(x1,y1,x2,y2)
	segments(x2-cap/2,y2,x2+cap/2,y2)
}

#Figure 6B
clades.sel<-c("Mammalia","Amniota")
colors<-c("#7453A2","#BAACD4","#F79520","#6EBE44","#C4E1AE")

for(sample in c("CS23_embryonic_cortex")) {

	pdf(paste("plot_aggregated_length_byCCE_",sample,"_selected.pdf",sep=""),width=10,height=5)
	par(las=2,cex=0.8,mfrow=c(1,2))
	
	for(num1 in c(3,4)) {
		clade1=clades[num1] #current clade
		if(clade1=="MtAmn") {
			data.sel<-data[grepl(sample,anno[,4]) & (grepl("Tetrapoda",anno[,5]) | grepl("Gnathostomata",anno[,5]) | grepl("Vertebrata",anno[,5]) ),]
		}
		else {
			data.sel<-data[grepl(sample,anno[,4]) & grepl(clade1,anno[,5]),]
		}
		
		plot(0,xlim=c(-2000,2000),ylim=c(0,7),pch="",main=clade1,axes=F,xlab="Aggregated Length",ylab="")
		
		for(num2 in 1:(num1+1)) {
			len.median<-median(data.sel[data.sel[,num2]>0,num2])
			len.quartiles<-quantile(data.sel[data.sel[,num2]>0,num2])
			rect(-len.median/2,num2-0.5,len.median/2,num2+0.5)
			#left
			#75 perc
			draw_horizontal_errbar(-len.median/2,num2,-(len.median/2+(len.quartiles[4]-len.median)/2),num2,0.2)
			#25 perc
			draw_horizontal_errbar(-len.median/2,num2,-(len.median/2-(len.median-len.quartiles[2])/2),num2,0.2)			
			#right
			#75 perc
			draw_horizontal_errbar(len.median/2,num2,len.median/2+(len.quartiles[4]-len.median)/2,num2,0.2)
			#25 perc
			draw_horizontal_errbar(len.median/2,num2,len.median/2-(len.median-len.quartiles[2])/2,num2,0.2)				
		}
		
		axis(1,at=seq(-2000,2000,500),labels=seq(-2000,2000,500))
		axis(2,at=1:6,labels=c("Whole Enhancer",clades))
		box()
		
	}
	dev.off()
	
}

#Figure 6A

for(sample in c("CS23_embryonic_cortex")) {

	#pdf(paste("plot_aggregated_length_byCCE_",sample,"_selected.pdf",sep=""),width=10,height=5)

	
	data.sel.len<-list()
	for(num1 in 1:length(clades)) {
		clade1=clades[num1] #current clade
		if(clade1=="MtAmn") {
			data.sel<-data[grepl(sample,anno[,4]) & (grepl("Tetrapoda",anno[,5]) | grepl("Gnathostomata",anno[,5]) | grepl("Vertebrata",anno[,5]) ),]
		}
		else {
			data.sel<-data[grepl(sample,anno[,4]) & grepl(clade1,anno[,5]),]
		}
		data.sel.len[[clade1]]<-data.sel[,1]
	}
	data.sel.len.quantiles<-lapply(data.sel.len,quantile)
	
	
	par(las=2,cex=0.8,mfrow=c(1,2))
		
	#boxplot without whiskers
	boxplot(data.sel.len,ylim=c(0,9000),outline=F, width=rep(0.3,5),col =colors)
	
	#dot plot with error bars
	
	plot(0,xlim=c(1,5), ylim=c(0,9000),ylab="Length",xlab="",axes=F,pch="")
	for(num in 1:5) {
		#75%
		draw_vertical_errbar(num,as.numeric(lapply(data.sel.len.quantiles,function(x){x[3]}))[num],num,as.numeric(lapply(data.sel.len.quantiles,function(x){x[4]})[num]),0.05)
		#25%
		draw_vertical_errbar(num,as.numeric(lapply(data.sel.len.quantiles,function(x){x[3]}))[num],num,as.numeric(lapply(data.sel.len.quantiles,function(x){x[2]})[num]),0.05)
	}
	points(cbind(1:5,lapply(data.sel.len.quantiles,function(x){x[3]})),pch=19,col=colors,cex=2)
	axis(1,at=1:5,labels=clades)
	axis(2,at=seq(0,8000,2000),labels=seq(0,8000,2000))
	box()
	
}
	
#violin
vioplot(data.sel.len[["Eutheria"]],data.sel.len[["Theria"]],data.sel.len[["Mammalia"]],data.sel.len[["Amniota"]],data.sel.len[["MtAmn"]])


#lapply(data.sel.len,quantile)
$Eutheria
   0%   25%   50%   75%  100% 
  500  1125  1700  2700 17425 

$Theria
   0%   25%   50%   75%  100% 
  550  1375  2175  3425 35375 

$Mammalia
   0%   25%   50%   75%  100% 
  525  1450  2300  3650 27750 

$Amniota
      0%      25%      50%      75%     100% 
  525.00  1406.25  2300.00  3825.00 32075.00 

$MtAmn
   0%   25%   50%   75%  100% 
  525  1675  2675  4500 39150 
  
#lapply(data.sel.len,length) (Correct, the same with Figure 1)
$Eutheria
[1] 2385

$Theria
[1] 2017

$Mammalia
[1] 2523

$Amniota
[1] 1770

$MtAmn
[1] 1395


