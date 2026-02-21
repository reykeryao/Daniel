rm(list=ls())
library(tidyverse)
library(RColorBrewer)
library(vioplot)
setwd("/stor/work/Lambowitz/yaojun/Work/Daniel/JA26026/MAFFT/")
'''
for (sample in paste0("DK6_",2:7)){
  print(paste0("Start processing Dataset ",sample))
  dat<-read.table(paste0(sample,".output"),
                  col.names=c("ID","Align"))
  dat$ID<-gsub("_","@",dat$ID)
  dat<-separate(dat,ID,into=c("ID","Reads"),sep="@",convert=T)
  dat<-separate(dat,Align,into=c("Empty",paste0("Pos",1:nchar(dat$Align[1]))),
                sep="")
  dat<-dat[,-3:-8]
  colnames(dat)[-1:-2]<-paste0("Pos",1:50)
  ### fill in position 1-5 and 46-50, these are not used in error rate and used to discard reads contain long deletions
  dat[-1,3:7]<-dat[1,3:7]
  dat[-1,48:52]<-dat[1,48:52]
  dat$Seq<-apply(dat[,-1:-2],1,function(x){unlist(paste(x,collapse = ""))})
  dat<-dat[!grepl("[A|C|G|T|N]-{2,}[A|C|G|T|N]",dat$Seq),]
  ## kepp all reads in all
  dat_all<-dat
  ### separate Pos31 based on A and T
  dat1<-dat[dat$Pos31=="A",]
  dat2<-dat[dat$Pos31!="A",]
  ### calculate the number of reads need to be removed
  ## if total corrected less than half of allreads, remove them all
  half_total<-round(sum(dat$Reads[-1])/2)
  total_correct<-sum(dat1$Reads[-1])
  if (total_correct<=half_total){
    dat<-rbind(dat[1,],dat2)
    dat_rm<-dat1
  } else {
    ## proportional reduced based on abundance
    dat_rm<-dat1
    dat1$Reads[-1]<-dat1$Reads[-1]*(1-(half_total/total_correct))
    dat<-rbind(dat1,dat2)
    dat_rm$Reads[-1]<-dat_rm$Reads[-1]*half_total/total_correct
  }
  ## check overall error rate of pos 6 to 45
  for (i in 8:47){
    if (i==8){
      tmp<-data.frame(table(rep(dat[-1,i],dat$Reads[-1])))
      tmp1<-data.frame(table(rep(dat_rm[-1,i],dat_rm$Reads[-1])))
      tmp_all<-tmp<-data.frame(table(rep(dat_all[-1,i],dat_all$Reads[-1])))
    } else {
      tmp<-merge(tmp, data.frame(table(rep(dat[-1,i],dat$Reads[-1]))),by=1,all=T)
      tmp1<-merge(tmp1, data.frame(table(rep(dat_rm[-1,i],dat_rm$Reads[-1]))),by=1,all=T)
      tmp_all<-merge(tmp_all, data.frame(table(rep(dat_all[-1,i],dat_all$Reads[-1]))),by=1,all=T)
    }
  }
  tmp_all[is.na(tmp_all)]<-0
  tmp[is.na(tmp)]<-0
  tmp1[is.na(tmp1)]<-0
  colnames(tmp)<-colnames(tmp1)<-colnames(tmp_all)<-c("NT",colnames(dat)[8:47])
  if (sample=="DK6_2"){
    res_all<-setNames(list(tmp_all),sample)
    res<-setNames(list(tmp),sample)
    res1<-setNames(list(tmp1),sample)
  } else {
    res_all<-c(res_all,setNames(list(tmp_all),sample))
    res<-c(res,setNames(list(tmp),sample))
    res1<-c(res1,setNames(list(tmp1),sample))
  }
}
names(res)
saveRDS(res_all,"processed_mafft.outout")
saveRDS(res,"processed_mafft_corrected.outout")
saveRDS(res1,"processed_mafft_removed.outout")
'''
res_name<-c("processed_mafft.outout","processed_mafft_corrected.outout","processed_mafft_removed.outout")

bcol<-c("lightblue","blue","gold","tomato","gray")
vcol<-brewer.pal(6,name = "Set1")
r_ran<-list(c(1:3,5:24,28:32),c(1:3,5:24),c(1:3,5:12),c(13:24),c(28:32))
names(r_ran)<-c("cDNA (all regions)","cDNA (region not coved by primer)",
                "cDNA (3' end)","cDNA (middle)","cDNA (5' end)")
for (sample in 1:length(res_name)){
  res<-readRDS(res_name[sample])
  for (i in 1:length(res)){
    tmp<-res[[i]]
    tmp_er<-apply(tmp[,-1],2,function(x){max(x)/sum(x)})
    if (i==1){
      Mismatch<-data.frame(tmp[,c("NT","Pos31")])
      Er<-data.frame(tmp_er)
    } else {
      Mismatch<-merge(Mismatch,data.frame(tmp[,c("NT","Pos31")]),by="NT",all=T)
      Er<-cbind(Er,data.frame(tmp_er))
    }
  }
  rownames(Mismatch)<-Mismatch[,1]
  rownames(Mismatch)[-1]<-c("T","G","C","N","A")
  Mismatch<-Mismatch[c("A","C","G","T","N"),-1]
  
  Mismatch<-Mismatch[rowSums(Mismatch)>0,]
  colnames(Mismatch)<-c("PPRT:Δ1184-1776","PPRT:Δ379-1776","GsI-IIC",
                        "GsI-IICxPPRT","HIV","HIVxPPRT")
  Mismatch<-prop.table(as.matrix(Mismatch),2)*100
  Mismatch_alt<-Mismatch
  Mismatch_alt["T",]<-Mismatch_alt["T",]-sapply(Mismatch_alt["T",],function(x){min(x,50)})
  Mismatch_alt<-prop.table(as.matrix(Mismatch_alt),2)*100
  
  colnames(Er)<-colnames(Mismatch)
  Er<-100-100*Er
  if (sample==2){
    Corrected_Er<-Er
  }
  ### make plots
  pdf(paste0("../Figs/JA26026_",sample,".pdf"),width=11,height=8)
  ### barplot
  par(mar=c(10,5,3,1),mfrow=c(2,3))
  if (sample==1){
    mp<-barplot(cbind(Mismatch,NA,c(50,0,0,50),c(0,0,0,100),NA,NA),names.arg = rep(NA,11),
                col=bcol,ylim=c(0,100),
                ylab="Reads",yaxt="n")
    axis(1,at=mp[8:9],labels= c("No correction","100% correction"),las=2)
  } else {
    mp<-barplot(cbind(Mismatch,NA,NA,NA,NA,NA),names.arg = rep(NA,11),
                col=bcol,ylim=c(0,100),
                ylab="Reads",yaxt="n")
  }
  legend("right",legend = rownames(Mismatch),fill=bcol,bty="n")
  axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
  axis(1,at=mp[1:6],labels= colnames(Mismatch),las=2)
  
  mtext(side=3,at=mean(mp[1:9]),line = 1,"Pos31 (A->T)")
  # plot each region error rate
  for (j in 1:length(r_ran)){
    tmp_plot<-Er[r_ran[[j]],]
    y_max<-0.5
    #ceiling(max(tmp_plot))
    plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,y_max),
         axes=FALSE,ann=FALSE)
    vioplot(tmp_plot,main=names(r_ran)[j],add=TRUE,col=vcol)
    axis(1,at=1:6,labels=colnames(tmp_plot),las=2)
    axis(2,at=c(0,y_max),labels=c(0,y_max),las=2)
    mtext(side=2,line = 2,"Error rate (%)")
    t1<-t.test(tmp_plot[,1],tmp_plot[,2])$p.value
    t2<-t.test(tmp_plot[,3],tmp_plot[,4])$p.value
    t3<-t.test(tmp_plot[,5],tmp_plot[,6])$p.value
    t1<-case_when(t1 < 0.001 ~ "***",t1 < 0.01 ~ "**",
                  t1 < 0.05 ~ "*",.default = "n.s.")
    t2<-case_when(t2 < 0.001 ~ "***",t2 < 0.01 ~ "**",
                  t2 < 0.05 ~ "*",.default = "n.s.")
    t3<-case_when(t3 < 0.001 ~ "***",t3 < 0.01 ~ "**",
                  t3 < 0.05 ~ "*",.default = "n.s.")
    segments(1,y_max*0.8,2,y_max*0.8)
    text(1.5,y_max*0.85,t1,cex=0.75)
    segments(3,y_max*0.8,4,y_max*0.8)
    text(3.5,y_max*0.85,t2,cex=0.75)
    segments(5,y_max*0.8,6,y_max*0.8)
    text(5.5,y_max*0.85,t3,cex=0.75)
    if(j==1){
      legend("topleft",legend = c("***: p<0.001","**: p<0.01",
                                  "*: p<0.05","n.s.: non-significant"),
             bty="n",cex=0.75)
    }
  }
  dev.off()
}
                              
### Error by preference (of template nucleotides)
dat<-read.table(paste0("DK6_2",".output"),col.names=c("ID","Align"))
dat<-separate(dat,ID,into=c("ID","Reads"),sep="@",convert=T)
dat<-separate(dat,Align,into=c("Empty",paste0("Pos",1:nchar(dat$Align[1]))),sep="")
template<-dat[1,14:53]
Er<-cbind(Corrected_Er,t(template))

## not covered error rate violin by nucleotides
pdf("../Figs/violin_by_nt.pdf",width=6,height=8)
par(mar=c(2,5,3,1),mfrow=c(3,2))
tmp<-Er[c(1:3,5:24),]
plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,1),
     axes=FALSE,ann=FALSE)
vioplot(list(tmp[tmp$`1`=="A",1],
          tmp[tmp$`1`=="C",1],
          tmp[tmp$`1`=="T",1],
        tmp[tmp$`1`=="A",2],
        tmp[tmp$`1`=="C",2],
        tmp[tmp$`1`=="T",2]),
        main="cDNA (region not coved by primer)",add=TRUE,col=vcol[1:3])
axis(1,at=1:6,labels=rep(c("A","C","T"),2),las=1)
axis(2,at=seq(0,1,1),labels=seq(0,1,1),las=2)
mtext(side=2,line = 2,"Error rate (%)")
segments(1,0.6,3,0.6)
text(2,0.7,"PPRT",cex=0.5)
segments(4,0.6,6,0.6)
text(5,0.7,"-Exo",cex=0.5)

plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=c(0,5),
     axes=FALSE,ann=FALSE)
vioplot(list(tmp[tmp$`1`=="A",2]/tmp[tmp$`1`=="A",1],
             tmp[tmp$`1`=="C",2]/tmp[tmp$`1`=="C",1],
             tmp[tmp$`1`=="T",2]/tmp[tmp$`1`=="T",1]),
        main="cDNA (region not coved by primer)",add=TRUE,col=vcol[1:3])
axis(1,at=1:3,labels=rep(c("A","C","T")),las=1)
axis(2,at=seq(0,5,1),labels=seq(0,5,1),las=2)
mtext(side=2,line = 2,"Fidelity (fold change)")


plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,1),
     axes=FALSE,ann=FALSE)
vioplot(list(tmp[tmp$`1`=="A",3],
             tmp[tmp$`1`=="C",3],
             tmp[tmp$`1`=="T",3],
             tmp[tmp$`1`=="A",4],
             tmp[tmp$`1`=="C",4],
             tmp[tmp$`1`=="T",4]),
        main="cDNA (region not coved by primer)",add=TRUE,col=vcol[1:3])
axis(1,at=1:6,labels=rep(c("A","C","T"),2),las=1)
axis(2,at=seq(0,1,1),labels=seq(0,1,1),las=2)
mtext(side=2,line = 2,"Error rate (%)")
segments(1,0.6,3,0.6)
text(2,0.7,"GsI-IIc",cex=0.5)
segments(4,0.6,6,0.6)
text(5,0.7,"GsI-IIC x PPRT",cex=0.5)

plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=c(0,5),
     axes=FALSE,ann=FALSE)
vioplot(list(tmp[tmp$`1`=="A",3]/tmp[tmp$`1`=="A",4],
             tmp[tmp$`1`=="C",3]/tmp[tmp$`1`=="C",4],
             tmp[tmp$`1`=="T",3]/tmp[tmp$`1`=="T",4]),
        main="cDNA (region not coved by primer)",add=TRUE,col=vcol[1:3])
axis(1,at=1:3,labels=rep(c("A","C","T")),las=1)
axis(2,at=seq(0,5,1),labels=seq(0,5,1),las=2)
mtext(side=2,line = 2,"Fidelity (fold change)")

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,1),
     axes=FALSE,ann=FALSE)
vioplot(list(tmp[tmp$`1`=="A",5],
             tmp[tmp$`1`=="C",5],
             tmp[tmp$`1`=="T",5],
             tmp[tmp$`1`=="A",6],
             tmp[tmp$`1`=="C",6],
             tmp[tmp$`1`=="T",6]),
        main="cDNA (region not coved by primer)",add=TRUE,col=vcol[1:3])
axis(1,at=1:6,labels=rep(c("A","C","T"),2),las=1)
axis(2,at=seq(0,1,1),labels=seq(0,1,1),las=2)
mtext(side=2,line = 2,"Error rate (%)")
segments(1,0.6,3,0.6)
text(2,0.7,"HIV RT",cex=0.5)
segments(4,0.6,6,0.6)
text(5,0.7,"HIV RT x PPRT",cex=0.5)

plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=c(0,5),
     axes=FALSE,ann=FALSE)
vioplot(list(tmp[tmp$`1`=="A",5]/tmp[tmp$`1`=="A",6],
             tmp[tmp$`1`=="C",5]/tmp[tmp$`1`=="C",6],
             tmp[tmp$`1`=="T",5]/tmp[tmp$`1`=="T",6]),
        main="cDNA (region not coved by primer)",add=TRUE,col=vcol[1:3])
axis(1,at=1:3,labels=rep(c("A","C","T")),las=1)
axis(2,at=seq(0,5,1),labels=seq(0,5,1),las=2)
mtext(side=2,line = 2,"Fidelity (fold change)")
dev.off()


pdf("../Figs/error_by_position.pdf",height=8,width=8)
par(mfrow=c(3,1))
plot(y=c(Er[c(1:24),1],NA,NA,NA,Er[c(28:32),1]),x=1:32,ylim=c(0,0.5),bty="n",xlim=c(1,32),pch=19,cex=0.5,axes=FALSE,xlab="cDNA positions",ylab="Error rate (%)")
points(y=c(Er[c(1:24),2],NA,NA,NA,Er[c(28:32),2]),x=1:32,col="red",pch=19,cex=0.5)
axis(1,at=seq(1,32,1),labels=FALSE)
axis(1,at=1:32,labels=Er$`1`[1:32],cex.axis=0.5)
abline(v=26)
axis(2,at=c(0,0.5),labels=c(0,0.5),cex.axis=0.5,las=2)
segments(1,0.45,12,0.45)
text(mean(1:12),0.48,labels = "3' end",cex=0.5)
segments(13,0.45,24,0.45)
text(mean(13:24),0.48,labels = "middle",cex=0.5)
segments(28,0.45,32,0.45)
text(mean(28:32),0.48,labels = "5' end",cex=0.5)
legend(30,0.45,legend = c("PPRT","-Exo"),pch=19,col=c("black","red"),bty="n",
       cex=0.5)

plot(y=c(Er[c(1:24),4],NA,NA,NA,Er[c(28:32),1]),x=1:32,ylim=c(0,0.5),bty="n",xlim=c(1,32),pch=19,cex=0.5,axes=FALSE,xlab="cDNA positions",ylab="Error rate (%)")
points(y=c(Er[c(1:24),3],NA,NA,NA,Er[c(28:32),2]),x=1:32,col="red",pch=19,cex=0.5)
axis(1,at=seq(1,32,1),labels=FALSE)
axis(1,at=1:32,labels=Er$`1`[1:32],cex.axis=0.5)
abline(v=26)
axis(2,at=c(0,0.5),labels=c(0,0.5),cex.axis=0.5,las=2)
segments(1,0.45,12,0.45)
text(mean(1:12),0.48,labels = "3' end",cex=0.5)
segments(13,0.45,24,0.45)
text(mean(13:23),0.48,labels = "middle",cex=0.5)
segments(28,0.45,32,0.45)
text(mean(28:32),0.48,labels = "5' end",cex=0.5)
legend(30,0.45,legend = c("GsI-IIC","GsI-IIC x PPRT"),pch=19,col=c("black","red"),bty="n",
       cex=0.5)

plot(y=c(Er[c(1:24),6],NA,NA,NA,Er[c(28:32),1]),x=1:32,ylim=c(0,0.5),bty="n",xlim=c(1,32),pch=19,cex=0.5,axes=FALSE,xlab="cDNA positions",ylab="Error rate (%)")
points(y=c(Er[c(1:24),5],NA,NA,NA,Er[c(28:32),2]),x=1:32,col="red",pch=19,cex=0.5)
axis(1,at=seq(1,32,1),labels=FALSE)
axis(1,at=1:32,labels=Er$`1`[1:32],cex.axis=0.5)
abline(v=26)
axis(2,at=c(0,0.5),labels=c(0,0.5),cex.axis=0.5,las=2)
segments(1,0.45,12,0.45)
text(mean(1:12),0.48,labels = "3' end",cex=0.5)
segments(13,0.45,24,0.45)
text(mean(13:23),0.48,labels = "middle",cex=0.5)
segments(28,0.45,32,0.45)
text(mean(28:32),0.48,labels = "5' end",cex=0.5)
legend(30,0.45,legend = c("HIV RT","HIV RT x PPRT"),pch=19,col=c("black","red"),bty="n",
       cex=0.5)
dev.off()


Er_mean<-rbind(colMeans(Er[28:32,-7]),colMeans(Er[13:23,-7]),colMeans(Er[1:12,-7]),
               colMeans(Er[1:23,-7]),colMeans(Er[c(1:23,28:32),-7]))
rownames(Er_mean)<-c("5' end","middle","3' end","not coved by primer","all regions")                                    


Er_mean

##combine 
res_name<-"processed_mafft_corrected.outout"

bcol<-c("lightblue","blue","gold","tomato","gray")
vcol<-brewer.pal(6,name = "Set1")
res<-readRDS(res_name)
for (i in 1:length(res)){
  tmp<-res[[i]]
  if (i==1){
    Mismatch<-data.frame(tmp[,c("NT","Pos31")])
  } else {
    Mismatch<-merge(Mismatch,data.frame(tmp[,c("NT","Pos31")]),by="NT",all=T)
  }
}

res<-readRDS(paste0("../../JA25159_211_284_330/JA25330/MAFFT/",res_name))
for (i in 1:length(res)){
  tmp<-res[[i]]
  Mismatch<-merge(Mismatch,data.frame(tmp[,c("NT","Pos31")]),by="NT",all=T)
}

Mismatch<-Mismatch[Mismatch$NT%in%c("A","C","G","T"),]
for (i in 2:7){
  Mismatch[,i]<-Mismatch[,i]+Mismatch[,i+6]
}
rownames(Mismatch)<-Mismatch$NT
rownames(Mismatch)<-c("T","G","C","A")
Mismatch<-prop.table(as.matrix(Mismatch[,2:7]),2)*100
colnames(Mismatch)<-c("PPRT:Δ1184-1776","PPRT:Δ379-1776","GsI-IIC",
                      "GsI-IICxPPRT","HIV","HIVxPPRT")
Mismatch<-Mismatch[c("A","C","G","T"),]
pdf("../Figs/Combined.pdf")
par(mar=c(8,5,3,3))
mp<-barplot(cbind(Mismatch,NA,NA),names.arg = rep(NA,8),
              col=bcol,ylim=c(0,100),
              ylab="Reads",yaxt="n")
axis(1,at=mp[1:6],labels= colnames(Mismatch),las=2)
legend("right",legend = rownames(Mismatch),fill=bcol,bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:6],labels= colnames(Mismatch),las=2)
dev.off()
