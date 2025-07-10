rm(list=ls())
library(tidyverse)
library(RColorBrewer)
setwd("/stor/work/Lambowitz/yaojun/Work/JA25159_25211_Daniel/Daniel/MAFFT/")
'''
### process MAFFT output
for (sample in c("DK1","DK2","DK3","DK4","DK5","DK6","DK7","DK8")){
  print(paste0("Start processing Dataset ",sample))
  dat<-read.table(paste0(sample,".output"),
                  col.names=c("ID","Align"))
  dat$ID<-gsub("_","@",dat$ID)
  dat<-separate(dat,ID,into=c("ID","Reads"),sep="@",convert=T)
  dat<-separate(dat,Align,into=c("Empty",paste0("Pos",1:nchar(dat$Align[1]))),
                sep="")
  dat<-dat[,-3]
  ### remove reads with too many truncations in 3 end (start of the primer)
  ### Pos50 can not be -
  dat<-dat[dat$Pos50!="-",]
  ### fill in the template sequence for Pos51 to Pos55, when calculate error rate
  ### Pos51 to Pos55 are not included
  dat$Pos51<-dat$Pos51[1]
  dat$Pos52<-dat$Pos52[1]
  dat$Pos53<-dat$Pos53[1]
  dat$Pos54<-dat$Pos54[1]
  dat$Pos55<-dat$Pos55[1]
  ### remake sequence
  dat$Seq<-apply(dat[,3:57],1,function(x){unlist(paste(x,collapse = ""))})
  dat$Seq<-gsub("-","",dat$Seq)
  dat$Len<-nchar(dat$Seq)
  print(dat$Seq[grep("N",dat$Seq)])
  if (sample=="DK1"){
    res<-setNames(list(dat),sample)
  } else {
    res<-c(res,setNames(list(dat),sample))
  }
}
saveRDS(res,"processed_mafft_JA25211.outout")
'''
res<-readRDS("processed_mafft_JA25211.outout")
for (i in 1:length(res)){
  tmp<-res[[i]]
  tmp<-tmp[tmp$Pos35!="-" & tmp$Pos36!="-",]
  tmp<-data.frame(table(rep(tmp$Len[-1],tmp$Reads[-1])))
  if (i==1){
    ttt<-merge(data.frame("Len"=c(10:60)),tmp,by=1,all=T)
  } else {
    ttt<-merge(ttt,tmp,by=1,all=T)
  }
}
colnames(ttt)[-1]<-names(res)
ttt[is.na(ttt)]<-0
ttt[,-1]<-100*prop.table(as.matrix(ttt[,-1]),2)
col=c("black","lightblue","tomato")
pdf("../Figs/cDNA_lenth_JA25211.pdf")
par(mfrow=c(3,1))
plot(ttt$DK1~ttt$Len,lwd=2,col=col[1],type="l",xlab="Length (nt)",ylim=c(0,100),
     main="cDNA length distribution")
lines(ttt$DK2~ttt$Len,lwd=2,col=col[2])
lines(ttt$DK3~ttt$Len,lwd=2,col=col[3])
legend("topleft",legend = c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA"),lty=1,lwd=2,
       bty="n",col=col)
abline(v=50,lty=2)
abline(v=19,lty=2)
text(50,80,"full length",pos=4)

plot(ttt$DK4~ttt$Len,lwd=2,type="l",xlab="Length (nt)",ylim=c(0,100),
     main="cDNA length distribution")
lines(ttt$DK5~ttt$Len,lwd=2,col="red")
legend("topleft",legend = c("GsI IIC","GIIxPPRT"),lty=1,lwd=2,
       bty="n",col=c("black","red"))
abline(v=50,lty=2)
abline(v=19,lty=2)
text(50,80,"full length",pos=4)

plot(ttt$DK6~ttt$Len,lwd=2,col=col[1],type="l",xlab="Length (nt)",ylim=c(0,100),
     main="cDNA length distribution")
lines(ttt$DK7~ttt$Len,lwd=2,col=col[2])
lines(ttt$DK8~ttt$Len,lwd=2,col=col[3])
legend("topleft",legend = c("HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),lty=1,lwd=2,
       bty="n",col=col)
abline(v=50,lty=2)
abline(v=19,lty=2)
text(50,80,"full length",pos=4)
dev.off()


### error rate and mismath at Pos36
for (i in 1:length(res)){
  tmp<-res[[i]]
  tmp<-tmp[tmp$Pos35!="-" & tmp$Pos36!="-",]
  tmp36<-tmp[tmp$Pos36!="-",]
  tmp36<-tmp36[-1,]
  err<-data.frame(apply(tmp[,8:57],2,function(x){
    temp<-x[1];
    pos<-rep(x[-1],tmp$Reads[-1])
    (sum(pos==temp)/sum(pos!="-"))
  }))
  if (i==1){
    res36<-merge(data.frame("NT"=c("A","C","G","T")),
                 data.frame(table(rep(tmp36$Pos36,tmp36$Reads))),by=1,all=T)
    eR<-err
  } else {
    res36<-merge(res36,
                 data.frame(table(rep(tmp36$Pos36,tmp36$Reads))),by=1,all=T)
    eR<-cbind(eR,err)
  }
}
names(res36)[-1]<-names(eR)<-names(res)
res36[is.na(res36)]<-0
res36[,-1]<-prop.table(as.matrix(res36[,-1]),2)*100
eR<-100-100*eR
### 33 and 36 pos bar
pdf("../Figs/mismatch_errorrate_JA25211.pdf")
par(mfrow=c(2,2))
frq<-res36
rownames(frq)<-frq$NT
frq<-as.matrix(frq[,-1])
rownames(frq)<-rev(rownames(frq))
frq<-frq[sort(rownames(frq)),]
mp<-barplot(cbind(frq,NA,NA),names.arg = rep(NA,10),
        col=c("lightblue","blue","gold","tomato"),ylim=c(0,100),
        ylab="Reads",main="Nucleotide frequency at Pos31 (A)",yaxt="n")
legend("right",legend = rownames(frq),fill=c("lightblue","blue","gold","tomato"),bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:8],labels= c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA",
                            "GsI IIC","GIIxPPRT","HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),las=2)

col=c("black","lightblue","tomato")
plot(eR$DK1[45:6]~c(5:44),type="l",ylim=c(0,100),xlab="cDNA (5'->3')",ylab="Error rate (%)")
lines(eR$DK2[45:6]~c(5:44),type="l",col=col[2])
lines(eR$DK3[45:6]~c(5:44),type="l",col=col[3])
legend("topleft",lty=1,col=col,bty="n",legend = c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA"))
abline(v=19,lty=2)
text(19,60,pos=4,"mismatched T->A")

plot(eR$DK4[45:6]~c(5:44),type="l",ylim=c(0,100),xlab="cDNA (5'->3')",ylab="Error rate (%)")
lines(eR$DK5[45:6]~c(5:44),type="l",col="red")
legend("topleft",legend = c("GsI IIC","GIIxPPRT"),lty=1,lwd=2,
       bty="n",col=c("black","red"))
abline(v=19,lty=2)
text(19,60,pos=4,"mismatched T->A")

plot(eR$DK6[45:6]~c(5:44),type="l",ylim=c(0,100),xlab="cDNA (5'->3')",ylab="Error rate (%)")
lines(eR$DK7[45:6]~c(5:44),type="l",col=col[2])
lines(eR$DK8[45:6]~c(5:44),type="l",col=col[3])
legend("topleft",legend = c("HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),lty=1,lwd=2,
       bty="n",col=col)
abline(v=19,lty=2)
text(19,60,pos=4,"mismatched T->A")
dev.off()

library(vioplot)
pdf("../Figs/error_rate_violinplot_JA25211.pdf",width=11,height=12)
vcol<-brewer.pal(8,"Paired")
par(mfrow=c(3,1))
vioplot(eR[32:45,],names=c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA",
                           "GsI IIC","GIIxPPRT","HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),
        main="cDNA (5' end)",ylab="Error rate (%)",
        ylim=c(0,2),col=vcol)
vioplot(eR[18:30,],names=c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA",
                                     "GsI IIC","GIIxPPRT","HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),
        main="cDNA (middle)",ylab="Error rate (%)",
        ylim=c(0,100),col=vcol)
vioplot(eR[6:17,],names=c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA",
                                    "GsI IIC","GIIxPPRT","HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),
        main="cDNA (3' end)",ylab="Error rate (%)",
        ylim=c(0,8),col=vcol)
dev.off()

### only FL

for (i in 1:length(res)){
  tmp<-res[[i]]
  tmp<-tmp[tmp$Pos11!="-",]
  err<-data.frame(apply(tmp[,8:57],2,function(x){
    temp<-x[1];
    pos<-rep(x[-1],tmp$Reads[-1])
    (sum(pos==temp)/sum(pos!="-"))
  }))
  if (i==1){
    eR<-err
  } else {
    eR<-cbind(eR,err)
  }
}
names(eR)<-names(res)
eR<-100-100*eR

pdf("../Figs/FL_error_rate_violinplot_JA25211.pdf",width=11,height=12)
vcol<-brewer.pal(8,"Paired")
par(mfrow=c(3,1))
vioplot(eR[32:45,],names=c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA",
                           "GsI IIC","GIIxPPRT","HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),
        main="cDNA (5' end)",ylab="Error rate (%)",
        ylim=c(0,8),col=vcol)
vioplot(eR[18:30,],names=c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA",
                           "GsI IIC","GIIxPPRT","HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),
        main="cDNA (middle)",ylab="Error rate (%)",
        ylim=c(0,8),col=vcol)
vioplot(eR[6:17,],names=c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA",
                          "GsI IIC","GIIxPPRT","HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),
        main="cDNA (3' end)",ylab="Error rate (%)",
        ylim=c(0,8),col=vcol)
dev.off()

### primer vs products
for (i in (1:length(res))){
  tmp<-res[[i]][-1,]
  tmp<-tmp[tmp$Pos36!="-",]
  pp<-data.frame("cDNA"=sum(tmp$Reads[tmp$Pos35!="-"]),
                 "Primer"=sum(tmp$Reads[tmp$Pos35=="-"]),
                 "Primer_A"=sum(tmp$Reads[tmp$Pos35=="-" & tmp$Pos36=="T"]),
                 "Primer_T"=sum(tmp$Reads[tmp$Pos35=="-" & tmp$Pos36=="A"]),
                 "cDNA_A"=sum(tmp$Reads[tmp$Pos35!="-" & tmp$Pos36=="T"]),
                 "cDNA_T"=sum(tmp$Reads[tmp$Pos35!="-" & tmp$Pos36=="A"]))
  if (i==1){
    ppt<-pp
  } else {
    ppt<-rbind(ppt,pp)
  }
}
rownames(ppt)<-names(res)
ppt<-t(ppt)

pdf("../Figs/Mismatched_JA25211.pdf",width=8,height=8)
par(mfrow=c(2,2),mar=c(8,3,2,1))
tmp<-ppt[1:2,]
tmp<-prop.table(tmp,2)*100
mp<-barplot(cbind(tmp,NA,NA,NA),names.arg = rep(NA,11),
            col=c("lightblue","tomato"),ylim=c(0,100),
            ylab="Reads",main="Extended primer (%)",yaxt="n")
legend("right",legend = c("Extended","Unused"),title = "Primer",fill=c("lightblue","tomato"),bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:8],labels= c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA",
                            "GsI IIC","GIIxPPRT","HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),las=2)

plot.new()
tmp<-ppt[3:4,]
tmp<-prop.table(tmp,2)*100
mp<-barplot(cbind(tmp,NA,NA,NA,NA),names.arg = rep(NA,12),
            col=c("lightblue","tomato"),ylim=c(0,100),
            ylab="Reads",main="Mismatched nulcotide in unused primer",yaxt="n")
legend("right",legend = c("Original (A)","Corrected (T)"),fill=c("lightblue","tomato"),bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:8],labels= c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA",
                             "GsI IIC","GIIxPPRT","HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),las=2)

tmp<-ppt[5:6,]
tmp<-prop.table(tmp,2)*100
mp<-barplot(cbind(tmp,NA,NA,NA,NA),names.arg = rep(NA,12),
            col=c("lightblue","tomato"),ylim=c(0,100),
            ylab="Reads",main="Mismatched nulcotide in extended primer",yaxt="n")
legend("right",legend = c("Original (A)","Corrected (T)"),fill=c("lightblue","tomato"),bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:8],labels= c("PPRT(D1184-1776)","PPRT(DD379-1776)","RT YVAA",
                            "GsI IIC","GIIxPPRT","HIV-1 RT","HIV-1xPPRT rep1","HIV-1xPPRT rep2"),las=2)
dev.off()
