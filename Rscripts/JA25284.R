rm(list=ls())
library(tidyverse)
library(RColorBrewer)
library(vioplot)
setwd("/stor/work/Lambowitz/yaojun/Work/JA25159_25211_Daniel/Daniel/JA25284/MAFFT/")
'''
### process MAFFT output
for (sample in paste0("DK3_",1:6)){
  print(paste0("Start processing Dataset ",sample))
  dat<-read.table(paste0(sample,".output"),
                  col.names=c("ID","Align"))
  dat$ID<-gsub("_","@",dat$ID)
  dat<-separate(dat,ID,into=c("ID","Reads"),sep="@",convert=T)
  dat<-separate(dat,Align,into=c("Empty",paste0("Pos",1:nchar(dat$Align[1]))),
                sep="")
  dat<-dat[,-3]
  ### fill in position 1-5 and 46-50, these are not used in error rate and used to discard reads contain long deletions
  dat[-1,3:7]<-dat[1,3:7]
  dat[-1,48:52]<-dat[1,48:52]
  dat$Seq<-apply(dat[,-1:-2],1,function(x){unlist(paste(x,collapse = ""))})
  dat<-dat[!grepl("[A|C|G|T|N]-{2,}[A|C|G|T|N]",dat$Seq),]
  ## check overall error rate of pos 6 to 45
  for (i in 8:47){
    if (i==8){
      tmp<-data.frame(table(rep(dat[-1,i],dat$Reads[-1])))
    } else {
      tmp<-merge(tmp, data.frame(table(rep(dat[-1,i],dat$Reads[-1]))),by=1,all=T)
    }
  }
  tmp[is.na(tmp)]<-0
  colnames(tmp)<-c("NT",colnames(dat)[8:47])
  if (sample=="DK3_1"){
    res<-setNames(list(tmp),sample)
  } else {
    res<-c(res,setNames(list(tmp),sample))
  }
}
names(res)
saveRDS(res,"processed_mafft.outout")
'''
res<-readRDS("processed_mafft.outout")
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
Mismatch<-Mismatch[c("A","C","G","T","N"),-1]

Mismatch<-Mismatch[rowSums(Mismatch)>0,]
colnames(Mismatch)<-c("PPRT:Δ1184-1776","PPRT:Δ379-1776","GsI-IIC","GsI-IICxPPRT","HIV","HIVxPPRT")
Mismatch<-prop.table(as.matrix(Mismatch),2)*100
Mismatch_alt<-Mismatch
Mismatch_alt["A",]<-Mismatch_alt["A",]-50
Mismatch_alt<-prop.table(as.matrix(Mismatch_alt),2)*100

colnames(Er)<-colnames(Mismatch)
Er<-100-100*Er
bcol<-c("lightblue","blue","gold","tomato","gray")
vcol<-brewer.pal(6,name = "Set1")


pdf("../Figs/JA25284.pdf",width=8,height=8)
### barplot
par(mar=c(10,5,3,1),mfcol=c(2,2))
mp<-barplot(cbind(Mismatch,NA),names.arg = rep(NA,7),
        col=bcol,ylim=c(0,100),
        ylab="Reads",yaxt="n")
legend("right",legend = rownames(Mismatch),fill=bcol,bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:6],labels= colnames(Mismatch),las=2)
mtext(side=3,at=mean(mp[1:6]),line = 1,"Pos31 (A->T)")

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,2),
     axes=FALSE,ann=FALSE)
vioplot(Er[28:40,],main="cDNA (5' end)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=c(0,1,2),labels=c(0,1,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){
#  for (j in (i+1):6){
#    print(t.test(Er[28:40,i],Er[28:40,j])$p.value)
#  }
#}
## all comb p value>0.1
plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,2),
     axes=FALSE,ann=FALSE)
vioplot(Er[13:24,],main="cDNA (middle)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=c(0,1,2),labels=c(0,1,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){
#    for (j in (i+1):6){
#      print(t.test(Er[13:24,i],Er[13:24,j])$p.value)
#    }
#  }
#[1] 0.3172193
#[1] 0.00149837
#[1] 0.0008919306
#[1] 0.0009951434
#[1] 0.0009168511
#[1] 2.941819e-05
#[1] 1.601766e-05
#[1] 1.943196e-05
#[1] 1.698258e-05
#[1] 0.07322565
#[1] 0.1460879
#[1] 0.08567111
#[1] 0.4468751
#[1] 0.8474687
#[1] 0.5527517

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,2),
     axes=FALSE,ann=FALSE)
vioplot(Er[1:12,],main="cDNA (3' end)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=c(0,1,2),labels=c(0,1,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")
#for (i in 1:5){
#    for (j in (i+1):6){
#      print(t.test(Er[13:24,i],Er[13:24,j])$p.value)
#    }
#  }
#[1] 0.3172193
#[1] 0.00149837
#[1] 0.0008919306
#[1] 0.0009951434
#[1] 0.0009168511
#[1] 2.941819e-05
#[1] 1.601766e-05
#[1] 1.943196e-05
#[1] 1.698258e-05
#[1] 0.07322565
#[1] 0.1460879
#[1] 0.08567111
#[1] 0.4468751
#[1] 0.8474687
#[1] 0.5527517

dev.off()

