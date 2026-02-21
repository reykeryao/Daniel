rm(list=ls())
library(tidyverse)
library(RColorBrewer)
library(vioplot)
setwd("/stor/work/Lambowitz/yaojun/Work/JA25159_25211_Daniel/JA25330/MAFFT/")
'''
### process MAFFT output
for (sample in paste0("DK4_",2:7)){
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
  if (sample=="DK4_2"){
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
rownames(Mismatch)[-1]<-c("T","G","C","N","A")
Mismatch<-Mismatch[c("A","C","G","T","N"),-1]

Mismatch<-Mismatch[rowSums(Mismatch)>0,]
colnames(Mismatch)<-c("PPRT:Δ1184-1776","PPRT:Δ379-1776","GsI-IIC","GsI-IICxPPRT","HIV","HIVxPPRT")
Mismatch<-prop.table(as.matrix(Mismatch),2)*100
Mismatch_alt<-Mismatch
Mismatch_alt["T",]<-Mismatch_alt["T",]-sapply(Mismatch_alt["T",],function(x){min(x,50)})
Mismatch_alt<-prop.table(as.matrix(Mismatch_alt),2)*100

colnames(Er)<-colnames(Mismatch)
Er<-100-100*Er
bcol<-c("lightblue","blue","gold","tomato","gray")
vcol<-brewer.pal(6,name = "Set1")


pdf("../Figs/JA25330.pdf",width=11,height=8)
### barplot
par(mar=c(10,5,3,1),mfrow=c(2,3))
mp<-barplot(cbind(Mismatch,NA,c(50,0,0,50),c(0,0,0,100),NA,NA),names.arg = rep(NA,11),
        col=bcol,ylim=c(0,100),
        ylab="Reads",yaxt="n")
legend("right",legend = rownames(Mismatch),fill=bcol,bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:6],labels= colnames(Mismatch),las=2)
axis(1,at=mp[8:9],labels= c("No correction","100% correction"),las=2)
mtext(side=3,at=mean(mp[1:9]),line = 1,"Pos31 (A->T)")
#
#mp<-barplot(cbind(Mismatch_alt,NA),names.arg = rep(NA,7),
#            col=bcol,ylim=c(0,100),
#            ylab="Reads",yaxt="n")
#legend("right",legend = rownames(Mismatch_alt),fill=bcol,bty="n")
#axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
#axis(1,at=mp[1:6],labels= colnames(Mismatch_alt),las=2)
#mtext(side=3,at=mean(mp[1:6]),line = 1,"Pos31 (A->T)")

plot.new()

# all region
plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,3),
     axes=FALSE,ann=FALSE)
vioplot(Er[c(1:24,28:40),],main="cDNA (all regions)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=c(0:3),labels=c(0:3),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){
#    for (j in (i+1):6){
#      print(t.test(Er[13:24,i],Er[13:24,j])$p.value)
#    }
#  }
#[1] 0.0006942058 *
#[1] 0.2688043
#[1] 0.07870656
#[1] 0.2250816
#[1] 0.06863153
#[1] 0.0005815936
#[1] 2.026358e-06
#[1] 3.182553e-05
#[1] 1.919138e-06
#[1] 0.005685414 *
#[1] 0.7562149
#[1] 0.004881747
#[1] 0.06939632
#[1] 0.2375967
#[1] 0.002283686 *

segments(1,2.9,2,2.9)
text(1.5,2.95,"***",cex=1)
segments(3,1.5,4,1.5)
text(3.5,1.55,"**",cex=1)
segments(5,0.8,6,0.8)
text(5.5,0.85,"**",cex=1)    
legend(4.5,3,legend = c("***: p<0.001","**: p<0.01","n.s.: non-significant"),bty="n",cex=0.75)

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,3),
     axes=FALSE,ann=FALSE)
vioplot(Er[1:12,],main="cDNA (3' end)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=c(0:3),labels=c(0:3),las=2)
mtext(side=2,line = 2,"Error rate (%)")
#for (i in 1:5){
#    for (j in (i+1):6){
#      print(t.test(Er[1:12,i],Er[1:12,j])$p.value)
#    }
#  }
#[1] 0.008747756
#[1] 0.5511241
#[1] 0.8330524
#[1] 0.2710094
#[1] 0.8680004
#[1] 0.01095252
#[1] 0.008002413
#[1] 0.01311614
#[1] 0.008166166
#[1] 0.4481272
#[1] 0.5853106
#[1] 0.4730048
#[1] 0.22748
#[1] 0.965639
#[1] 0.2414688
segments(1,2.9,2,2.9)
text(1.5,2.95,"***")
segments(3,0.8,4,0.8)
text(3.5,0.85,"n.s.",cex=0.75)
segments(5,0.8,6,0.8)
text(5.5,0.85,"n.s.",cex=0.75)


plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,3),
     axes=FALSE,ann=FALSE)
vioplot(Er[13:23,],main="cDNA (middle)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=c(0:3),labels=c(0:3),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){
#    for (j in (i+1):6){
#      print(t.test(Er[13:24,i],Er[13:24,j])$p.value)
#    }
#  }
#[1] 0.0005496795
#[1] 0.1326321
#[1] 0.2713674
#[1] 0.01246073
#[1] 0.2370567
#[1] 0.001478164
#[1] 0.0004425972
#[1] 0.0009932534
#[1] 0.0004360198
#[1] 0.07049558
#[1] 0.5433648
#[1] 0.06724209
#[1] 0.0009934376
#[1] 0.9352771
#[1] 0.0007989979

segments(1,2.9,2,2.9)
text(1.5,2.95,"***")
segments(3,1.5,4,1.5)
text(3.5,1.55,"p=0.07",cex=0.75)
segments(5,0.8,6,0.8)
text(5.5,0.85,"***")

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,3),
     axes=FALSE,ann=FALSE)
vioplot(Er[28:40,],main="cDNA (5' end)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=c(0:3),labels=c(0:3),las=2)
mtext(side=2,line = 2,"Error rate (%)")
segments(1,1,6,1)
text(3.5,1.1,"n.s.",cex=0.75)
#for (i in 1:5){
#  for (j in (i+1):6){
#    print(t.test(Er[28:40,i],Er[28:40,j])$p.value)
#  }
#}
## all comb p value>0.1
dev.off()

Er_mean<-rbind(colMeans(Er[28:40,]),colMeans(Er[13:24,]),colMeans(Er[1:12,]),
               colMeans(Er[1:24,]),colMeans(Er[c(1:24,28:40),]))
rownames(Er_mean)<-c("5' end","middle","3' end","not coved by primer","all regions")                                    


#middle and 3'
plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,3),
     axes=FALSE,ann=FALSE)
vioplot(Er[1:24,],main="cDNA (region not coved by primer)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=c(0:3),labels=c(0:3),las=2)
mtext(side=2,line = 2,"Error rate (%)")
#for (i in 1:5){
#    for (j in (i+1):6){
#      print(t.test(Er[1:12,i],Er[1:12,j])$p.value)
#    }
#  }
#[1] 0.008876004 *
#[1] 0.4906076
#[1] 0.7512679
#[1] 0.1678343
#[1] 0.7903798
#[1] 1.643449e-05
#[1] 8.714869e-06
#[1] 2.399802e-05
#[1] 9.074542e-06
#[1] 0.1391561 *
#[1] 0.6116875
#[1] 0.1474423
#[1] 0.06901328
#[1] 0.8891532
#[1] 0.1243832 *
segments(1,2.95,2,2.95)
text(1.5,3,"***")
segments(3,1.4,4,1.4)
text(3.5,1.5,"n.s.",cex=0.75)
segments(5,0.8,6,0.8)
text(5.5,0.9,"n.s.",cex=0.75)

plot(y=Er[,1],x=51-(6:45),ylim=c(0,3),bty="n",xlim=c(0,50),pch=19,cex=0.5,axes=FALSE,
     xlab="cDNA positions",ylab="Error rate (%)")
points(y=Er[,2],x=51-(6:45),col="red",pch=19,cex=0.5)
axis(1,at=seq(0,50,1),labels=FALSE)
axis(1,at=seq(0,50,5),labels=seq(0,50,5),cex.axis=0.5)
abline(v=20)


'''
### try random remove up to 50% corrected pos31 reads
### process MAFFT output
for (sample in paste0("DK4_",2:7)){
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
    } else {
      tmp<-merge(tmp, data.frame(table(rep(dat[-1,i],dat$Reads[-1]))),by=1,all=T)
      tmp1<-merge(tmp1, data.frame(table(rep(dat_rm[-1,i],dat_rm$Reads[-1]))),by=1,all=T)
    }
  }
  tmp[is.na(tmp)]<-0
  tmp1[is.na(tmp1)]<-0
  colnames(tmp)<-colnames(tmp1)<-c("NT",colnames(dat)[8:47])
  if (sample=="DK4_2"){
    res<-setNames(list(tmp),sample)
    res1<-setNames(list(tmp1),sample)
  } else {
    res<-c(res,setNames(list(tmp),sample))
    res1<-c(res1,setNames(list(tmp1),sample))
  }
}
names(res)
saveRDS(res,"processed_mafft_corrected.outout")
saveRDS(res1,"processed_mafft_removed.outout")
'''

res<-readRDS("processed_mafft_corrected.outout")
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
Mismatch<-Mismatch[c(-1,-5),]
rownames(Mismatch)<-c("T","G","C","A")
Mismatch<-Mismatch[c("A","C","G","T"),-1]

Mismatch<-Mismatch[rowSums(Mismatch)>0,]
colnames(Mismatch)<-c("PPRT:Δ1184-1776","PPRT:Δ379-1776","GsI-IIC","GsI-IICxPPRT","HIV","HIVxPPRT")
Mismatch<-prop.table(as.matrix(Mismatch),2)*100

colnames(Er)<-colnames(Mismatch)
Er<-100-100*Er
bcol<-c("lightblue","blue","gold","tomato","gray")
vcol<-brewer.pal(6,name = "Set1")

pdf("../Figs/JA25330_2.pdf",width=11,height=8)
### barplot
par(mar=c(10,5,3,1),mfrow=c(2,3))
mp<-barplot(cbind(Mismatch,NA),names.arg = rep(NA,7),
            col=bcol,ylim=c(0,100),
            ylab="Reads",yaxt="n")
legend("right",legend = rownames(Mismatch),fill=bcol,bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:6],labels= colnames(Mismatch),las=2)
mtext(side=3,at=mean(mp[1:9]),line = 1,"Pos31 (A->T)")

# all region
plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(Er[c(1:23,28:40),],main="cDNA (all regions)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){for (j in (i+1):6){print(t.test(Er[c(1:23,28:40),i],Er[c(1:23,28:40),j])$p.value)}}
#    
#[1] 2.879374e-05  ***
#[1] 0.006545534 **
#[1] 0.008286348 **

segments(1,5.5,2,5.5)
text(1.5,5.75,"***",cex=1)
segments(3,1.6,4,1.6)
text(3.5,1.85,"**",cex=1)
segments(5,1.6,6,1.6)
text(5.5,1.85,"**",cex=1)    
legend(4.5,5.8,legend = c("***: p<0.001","**: p<0.01","*: p<0.05","n.s.: non-significant"),bty="n",cex=0.75)

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(Er[1:12,],main="cDNA (3' end)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")
#for (i in 1:5){for (j in (i+1):6){print(t.test(Er[c(1:12),i],Er[c(1:12),j])$p.value)}}
#[1] 0.006211696 **
#[1] 0.08031285
#[1] 0.03176017

segments(1,4.3,2,4.3)
text(1.5,4.55,"**")
segments(3,1.2,4,1.2)
text(3.5,1.45,"n.s.",cex=1)
segments(5,1.2,6,1.2)
text(5.5,1.45,"*",cex=1)   


plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(Er[13:23,],main="cDNA (middle)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){for (j in (i+1):6){print(t.test(Er[c(13:23),i],Er[c(13:23),j])$p.value)}}
#[1] 0.0009942524   ***
#[1] 8.729678e-06  ***
#[1] 4.806121e-05 ***

segments(1,5.4,2,5.4)
text(1.5,5.65,"***")
segments(3,1.2,4,1.2)
text(3.5,1.45,"***",cex=1)
segments(5,1.2,6,1.2)
text(5.5,1.45,"***",cex=1)   

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(Er[28:40,],main="cDNA (5' end)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){for (j in (i+1):6){print(t.test(Er[c(28:40),i],Er[c(28:40),j])$p.value)}}
## 0.4755744  n.s
#[1] 0.009080529
#[1] 0.008194513  **
#[1] 0.01244064 *
segments(1,1.4,6,1.4)
text(3.5,1.65,"n.s.",cex=0.75)

#middle and 3'
plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(Er[1:23,],main="cDNA (region not coved by primer)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){for (j in (i+1):6){print(t.test(Er[c(1:23),i],Er[c(1:23),j])$p.value)}}
#[1] 1.418418e-05 ***
#[1] 0.0001248157  ***
#[1] 3.263295e-05 ***

segments(1,5.4,2,5.4)
text(1.5,5.65,"***")
segments(3,1.2,4,1.2)
text(3.5,1.45,"***",cex=1)
segments(5,1.2,6,1.2)
text(5.5,1.45,"***",cex=1)   

dev.off()

Er_mean<-rbind(colMeans(Er[28:40,]),colMeans(Er[13:23,]),colMeans(Er[1:12,]),
               colMeans(Er[1:23,]),colMeans(Er[c(1:23,28:40),]))
rownames(Er_mean)<-c("5' end","middle","3' end","not coved by primer","all regions")                                    
### Error by preference (of template nucleotides)
dat<-read.table(paste0("DK4_7",".output"),col.names=c("ID","Align"))
dat<-separate(dat,ID,into=c("ID","Reads"),sep="@",convert=T)
dat<-separate(dat,Align,into=c("Empty",paste0("Pos",1:nchar(dat$Align[1]))),sep="")
template<-dat[1,9:48]
Er<-cbind(Er,t(template))

## not covered error rate violin by nucleotides
pdf("../Figs/violin_by_nt.pdf",width=6,height=8)
par(mar=c(2,5,3,1),mfrow=c(3,2))
tmp<-Er[c(1:23),]
plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(list(tmp[tmp$`1`=="A",1],
          tmp[tmp$`1`=="C",1],
          tmp[tmp$`1`=="T",1],
        tmp[tmp$`1`=="A",2],
        tmp[tmp$`1`=="C",2],
        tmp[tmp$`1`=="T",2]),
        main="cDNA (region not coved by primer)",add=TRUE,col=vcol[1:3])
axis(1,at=1:6,labels=rep(c("A","C","T"),2),las=1)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")
segments(1,5.4,3,5.4)
text(2,5.6,"PPRT",cex=0.5)
segments(4,5.4,6,5.4)
text(5,5.6,"-Exo",cex=0.5)

plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=c(0,30),
     axes=FALSE,ann=FALSE)
vioplot(list(tmp[tmp$`1`=="A",2]/tmp[tmp$`1`=="A",1],
             tmp[tmp$`1`=="C",2]/tmp[tmp$`1`=="C",1],
             tmp[tmp$`1`=="T",2]/tmp[tmp$`1`=="T",1]),
        main="cDNA (region not coved by primer)",add=TRUE,col=vcol[1:3])
axis(1,at=1:3,labels=rep(c("A","C","T")),las=1)
axis(2,at=seq(0,30,5),labels=seq(0,30,5),las=2)
mtext(side=2,line = 2,"Fidelity (fold change")


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
axis(2,at=seq(0,1,0.5),labels=seq(0,1,0.5),las=2)
mtext(side=2,line = 2,"Error rate (%)")
segments(1,0.7,3,0.7)
text(2,0.8,"GsI-IIc",cex=0.5)
segments(4,0.7,6,0.7)
text(5,0.8,"GsI-IIC x PPRT",cex=0.5)

plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=c(0,10),
     axes=FALSE,ann=FALSE)
vioplot(list(tmp[tmp$`1`=="A",3]/tmp[tmp$`1`=="A",4],
             tmp[tmp$`1`=="C",3]/tmp[tmp$`1`=="C",4],
             tmp[tmp$`1`=="T",3]/tmp[tmp$`1`=="T",4]),
        main="cDNA (region not coved by primer)",add=TRUE,col=vcol[1:3])
axis(1,at=1:3,labels=rep(c("A","C","T")),las=1)
axis(2,at=seq(0,10,2),labels=seq(0,10,2),las=2)
mtext(side=2,line = 2,"Fidelity (fold change")

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
axis(2,at=seq(0,1,0.5),labels=seq(0,1,0.5),las=2)
mtext(side=2,line = 2,"Error rate (%)")
segments(1,0.7,3,0.7)
text(2,0.8,"HIV RT",cex=0.5)
segments(4,0.8,6,0.8)
text(5,0.9,"HIV RT x PPRT",cex=0.5)

plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=c(0,10),
     axes=FALSE,ann=FALSE)
vioplot(list(tmp[tmp$`1`=="A",5]/tmp[tmp$`1`=="A",6],
             tmp[tmp$`1`=="C",5]/tmp[tmp$`1`=="C",6],
             tmp[tmp$`1`=="T",5]/tmp[tmp$`1`=="T",6]),
        main="cDNA (region not coved by primer)",add=TRUE,col=vcol[1:3])
axis(1,at=1:3,labels=rep(c("A","C","T")),las=1)
axis(2,at=seq(0,10,2),labels=seq(0,10,2),las=2)
mtext(side=2,line = 2,"Fidelity (fold change")

dev.off()


pdf("../Figs/error_by_position.pdf",height=8,width=8)
par(mfrow=c(3,1))
plot(y=c(Er[c(1:23),1],NA,NA,NA,NA,Er[c(28:40),1]),x=1:40,ylim=c(0,6),bty="n",xlim=c(1,40),pch=19,cex=0.5,axes=FALSE,xlab="cDNA positions",ylab="Error rate (%)")
points(y=c(Er[c(1:23),2],NA,NA,NA,NA,Er[c(28:40),2]),x=1:40,col="red",pch=19,cex=0.5)
axis(1,at=seq(1,40,1),labels=FALSE)
axis(1,at=1:40,labels=Er$`1`,cex.axis=0.5)
abline(v=26)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),cex.axis=0.5,las=2)
segments(1,5.5,12,5.5)
text(mean(1:12),5.8,labels = "3' end",cex=0.5)
segments(13,5.5,23,5.5)
text(mean(13:23),5.8,labels = "middle",cex=0.5)
segments(28,5.5,40,5.5)
text(mean(28:40),5.8,labels = "5' end",cex=0.5)
legend(30,5,legend = c("PPRT","-Exo"),pch=19,col=c("black","red"),bty="n",
       cex=0.5)

plot(y=c(Er[c(1:23),4],NA,NA,NA,NA,Er[c(28:40),1]),x=1:40,ylim=c(0,6),bty="n",xlim=c(1,40),pch=19,cex=0.5,axes=FALSE,xlab="cDNA positions",ylab="Error rate (%)")
points(y=c(Er[c(1:23),3],NA,NA,NA,NA,Er[c(28:40),2]),x=1:40,col="red",pch=19,cex=0.5)
axis(1,at=seq(1,40,1),labels=FALSE)
axis(1,at=1:40,labels=Er$`1`,cex.axis=0.5)
abline(v=26)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),cex.axis=0.5,las=2)
segments(1,5.5,12,5.5)
text(mean(1:12),5.8,labels = "3' end",cex=0.5)
segments(13,5.5,23,5.5)
text(mean(13:23),5.8,labels = "middle",cex=0.5)
segments(28,5.5,40,5.5)
text(mean(28:40),5.8,labels = "5' end",cex=0.5)
legend(30,5,legend = c("GsI-IIC","GsI-IIC x PPRT"),pch=19,col=c("black","red"),bty="n",
       cex=0.5)

plot(y=c(Er[c(1:23),6],NA,NA,NA,NA,Er[c(28:40),1]),x=1:40,ylim=c(0,6),bty="n",xlim=c(1,40),pch=19,cex=0.5,axes=FALSE,xlab="cDNA positions",ylab="Error rate (%)")
points(y=c(Er[c(1:23),5],NA,NA,NA,NA,Er[c(28:40),2]),x=1:40,col="red",pch=19,cex=0.5)
axis(1,at=seq(1,40,1),labels=FALSE)
axis(1,at=1:40,labels=Er$`1`,cex.axis=0.5)
abline(v=26)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),cex.axis=0.5,las=2)
segments(1,5.5,12,5.5)
text(mean(1:12),5.8,labels = "3' end",cex=0.5)
segments(13,5.5,23,5.5)
text(mean(13:23),5.8,labels = "middle",cex=0.5)
segments(28,5.5,40,5.5)
text(mean(28:40),5.8,labels = "5' end",cex=0.5)
legend(30,5,legend = c("HIV RT","HIV RT x PPRT"),pch=19,col=c("black","red"),bty="n",
       cex=0.5)
dev.off()

### removed part
res<-readRDS("processed_mafft_removed.outout")
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
Mismatch<-Mismatch[c(-1,-5),]
rownames(Mismatch)<-c("T","G","C","A")
Mismatch<-Mismatch[c("A","C","G","T"),-1]

Mismatch<-Mismatch[rowSums(Mismatch)>0,]
colnames(Mismatch)<-c("PPRT:Δ1184-1776","PPRT:Δ379-1776","GsI-IIC","GsI-IICxPPRT","HIV","HIVxPPRT")
Mismatch<-prop.table(as.matrix(Mismatch),2)*100

colnames(Er)<-colnames(Mismatch)
Er<-100-100*Er
bcol<-c("lightblue","blue","gold","tomato","gray")
vcol<-brewer.pal(6,name = "Set1")

pdf("../Figs/JA25330_3.pdf",width=11,height=8)
### barplot
par(mar=c(10,5,3,1),mfrow=c(2,3))
mp<-barplot(cbind(Mismatch,NA),names.arg = rep(NA,7),
            col=bcol,ylim=c(0,100),
            ylab="Reads",yaxt="n")
legend("right",legend = rownames(Mismatch),fill=bcol,bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:6],labels= colnames(Mismatch),las=2)
mtext(side=3,at=mean(mp[1:9]),line = 1,"Pos31 (A->T)")

# all region
plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(Er[c(1:23,28:40),],main="cDNA (all regions)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){for (j in (i+1):6){print(t.test(Er[c(1:23,28:40),i],Er[c(1:23,28:40),j])$p.value)}}
#    
#[1] 0.004441044  **
#[1] 0.9667615 n.s.
#[1] 0.5575018 n.s.

segments(1,1.8,2,1.8)
text(1.5,2.05,"**",cex=1)
segments(3,1.8,4,1.8)
text(3.5,2.05,"n.s.",cex=1)
segments(5,1.8,6,1.8)
text(5.5,2.05,"n.s.",cex=1)    
legend(0.5,5.8,legend = c("**: p<0.01","*: p<0.05","n.s.: non-significant"),bty="n",cex=0.75)

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(Er[1:12,],main="cDNA (3' end)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")
#for (i in 1:5){for (j in (i+1):6){print(t.test(Er[c(1:12),i],Er[c(1:12),j])$p.value)}}
#[1] 0.03702697 *
#[1] 0.8182884
#[1] 0.9365555

segments(1,1.8,2,1.8)
text(1.5,2.05,"*",cex=1)
segments(3,1.8,4,1.8)
text(3.5,2.05,"n.s.",cex=1)
segments(5,1.8,6,1.8)
text(5.5,2.05,"n.s.",cex=1)    


plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(Er[13:23,],main="cDNA (middle)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){for (j in (i+1):6){print(t.test(Er[c(13:23),i],Er[c(13:23),j])$p.value)}}
#[1] 0.008596447   **
#[1] 0.7077861 
#[1] 0.8335082

segments(1,1.8,2,1.8)
text(1.5,2.05,"**",cex=1)
segments(3,1.8,4,1.8)
text(3.5,2.05,"n.s.",cex=1)
segments(5,1.8,6,1.8)
text(5.5,2.05,"n.s.",cex=1)   

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(Er[28:40,],main="cDNA (5' end)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){for (j in (i+1):6){print(t.test(Er[c(28:40),i],Er[c(28:40),j])$p.value)}}
## 0.4344592  n.s
#[1] 0.8772862  
#[1] 0.4829591
segments(1,1.8,2,1.8)
text(1.5,2.05,"n.s.",cex=1)
segments(3,1.8,4,1.8)
text(3.5,2.05,"n.s.",cex=1)
segments(5,1.8,6,1.8)
text(5.5,2.05,"n.s.",cex=1)   

#middle and 3'
plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,6),
     axes=FALSE,ann=FALSE)
vioplot(Er[1:23,],main="cDNA (region not coved by primer)",add=TRUE,col=vcol)
axis(1,at=1:6,labels=colnames(Er),las=2)
axis(2,at=seq(0,6,2),labels=seq(0,6,2),las=2)
mtext(side=2,line = 2,"Error rate (%)")

#for (i in 1:5){for (j in (i+1):6){print(t.test(Er[c(1:23),i],Er[c(1:23),j])$p.value)}}
#[1] 0.004213069 ***
#[1] 0.7232456
#[1] 0.9984317

segments(1,1.8,2,1.8)
text(1.5,2.05,"**",cex=1)
segments(3,1.8,4,1.8)
text(3.5,2.05,"n.s.",cex=1)
segments(5,1.8,6,1.8)
text(5.5,2.05,"n.s.",cex=1)  

dev.off()

Er_mean<-rbind(colMeans(Er[28:40,]),colMeans(Er[13:23,]),colMeans(Er[1:12,]),
               colMeans(Er[1:23,]),colMeans(Er[c(1:23,28:40),]))
rownames(Er_mean)<-c("5' end","middle","3' end","not coved by primer","all regions")                                    



plot(y=Er[,1],x=51-(6:45),ylim=c(0,6),bty="n",xlim=c(0,50),pch=19,cex=0.5,axes=FALSE,
     xlab="cDNA positions",ylab="Error rate (%)")
points(y=Er[,2],x=51-(6:45),col="red",pch=19,cex=0.5)
axis(1,at=seq(0,50,1),labels=FALSE)
axis(1,at=seq(0,50,5),labels=seq(0,50,5),cex.axis=0.5)
abline(v=20)
segments(5,5.5,17,5.5)
segments(22,5.5,32,5.5)
segments(33,5.5,45,5.5)





