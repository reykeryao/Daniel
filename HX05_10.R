rm(list=ls())
library(tidyverse)
library(RColorBrewer)
setwd("/stor/work/Lambowitz/yaojun/Work/JA25159_Daniel/MAFFT/")
'''
### process MAFFT output
for (sample in c("HX05","HX06","HX07","HX08","HX09","HX10")){
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
  if (sample=="HX05"){
    res<-setNames(list(dat),sample)
  } else {
    res<-c(res,setNames(list(dat),sample))
  }
}
saveRDS(res,"processed_mafft.outout")
'''
res<-readRDS("processed_mafft.outout")
for (i in 1:length(res)){
  tmp<-res[[i]]
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
col=c("black","lightblue","gold","tomato")
pdf("../Figs/cDNA_lenth.pdf")
par(mfrow=c(2,1))
plot(ttt$HX05~ttt$Len,lwd=2,col=col[1],type="l",xlab="Length (nt)",ylim=c(0,100),
     main="cDNA length distribution")
lines(ttt$HX06~ttt$Len,lwd=2,col=col[2])
lines(ttt$HX09~ttt$Len,lwd=2,col=col[3])
lines(ttt$HX10~ttt$Len,lwd=2,col=col[4])
legend("topleft",legend = c("PPRT","PrimPol","HIV RT","HIV-PPRT"),lty=1,lwd=2,
       bty="n",col=col)
abline(v=50,lty=2)
abline(v=19,lty=2)


plot(ttt$HX07~ttt$Len,lwd=2,type="l",xlab="Length (nt)",ylim=c(0,100),
     main="cDNA length distribution")
lines(ttt$HX08~ttt$Len,lwd=2,col="red")
legend("topleft",legend = c("PPRT","PrimPol"),lty=1,lwd=2,
       bty="n",col=c("black","red"))
abline(v=50,lty=2)
abline(v=23,lty=2)
text(23,40,"8-oxo-G",pos=2)
text(50,80,"full length",pos=4)
dev.off()


### error rate and mismath at Pos33 and Pos36
for (i in 1:length(res)){
  tmp<-res[[i]]
  tmp33<-tmp[tmp$Pos33!="-",]
  tmp33<-tmp33[-1,]
  tmp36<-tmp[tmp$Pos36!="-",]
  tmp36<-tmp36[-1,]
  err<-data.frame(apply(tmp[,8:57],2,function(x){
    temp<-x[1];
    pos<-rep(x[-1],tmp$Reads[-1])
    (sum(pos==temp)/sum(pos!="-"))
  }))
  if (i==1){
    res33<-merge(data.frame("NT"=c("A","C","G","T")),
                 data.frame(table(rep(tmp33$Pos33,tmp33$Reads))),by=1,all=T)
    res36<-merge(data.frame("NT"=c("A","C","G","T")),
                 data.frame(table(rep(tmp36$Pos36,tmp36$Reads))),by=1,all=T)
    eR<-err
  } else {
    res33<-merge(res33,
                 data.frame(table(rep(tmp33$Pos33,tmp33$Reads))),by=1,all=T)
    res36<-merge(res36,
                 data.frame(table(rep(tmp36$Pos36,tmp36$Reads))),by=1,all=T)
    eR<-cbind(eR,err)
  }
}
names(res33)[-1]<-names(res36)[-1]<-names(eR)<-names(res)
res33[is.na(res33)]<-0
res36[is.na(res36)]<-0
res33[,-1]<-prop.table(as.matrix(res33[,-1]),2)*100
res36[,-1]<-prop.table(as.matrix(res36[,-1]),2)*100
eR<-100-100*eR
### 33 and 36 pos bar
pdf("../Figs/mismatch_errorrate.pdf")
par(mfrow=c(3,2))
frq<-res33
rownames(frq)<-frq$NT
frq<-as.matrix(frq[,-1])
##RC rownames
rownames(frq)<-rev(rownames(frq))
frq<-frq[sort(rownames(frq)),]
mp<-barplot(cbind(frq,NA,NA),names.arg = rep(NA,8),
        col=c("lightblue","blue","gold","tomato"),ylim=c(0,100),
        ylab="Reads",main="Nucleotide frequency at Pos28 (8-oxo-G)",yaxt="n")
legend("right",legend = rownames(frq),fill=c("lightblue","blue","gold","tomato"),bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:6],labels= c("PPRT","Mut","PPRT|8-oxo-G",
              "Mut|8-oxo-G","HIV RT","HIV PPRT"),las=2)

frq<-res36
rownames(frq)<-frq$NT
frq<-as.matrix(frq[,-1])
rownames(frq)<-rev(rownames(frq))
frq<-frq[sort(rownames(frq)),]
mp<-barplot(cbind(frq,NA,NA),names.arg = rep(NA,8),
        col=c("lightblue","blue","gold","tomato"),ylim=c(0,100),
        ylab="Reads",main="Nucleotide frequency at Pos31 (A)",yaxt="n")
legend("right",legend = rownames(frq),fill=c("lightblue","blue","gold","tomato"),bty="n")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
axis(1,at=mp[1:6],labels= c("PPRT","Mut","PPRT|8-oxo-G",
                            "Mut|8-oxo-G","HIV RT","HIV PPRT"),las=2)

plot(eR$HX05[45:6]~c(5:44),type="l",ylim=c(0,100),xlab="cDNA (5'->3')",ylab="Error rate (%)")
lines(eR$HX06[45:6]~c(5:44),type="l",col="red")
legend("topleft",lty=1,col=c("black","red"),bty="n",legend = c("PPRT","PrimPol"))
abline(v=19,lty=2)
text(19,60,pos=4,"mismatched T->A")

plot(eR$HX09[45:6]~c(5:44),type="l",ylim=c(0,100),xlab="cDNA (5'->3')",ylab="Error rate (%)")
lines(eR$HX10[45:6]~c(5:44),type="l",col="red")
legend("topleft",lty=1,col=c("black","red"),bty="n",legend = c("HIV RT","HIV PPRT"))
abline(v=19,lty=2)
text(19,60,pos=4,"mismatched T->A")

plot(eR$HX07[45:6]~c(5:44),type="l",ylim=c(0,100),xlab="cDNA (5'->3')",ylab="Error rate (%)")
lines(eR$HX08[45:6]~c(5:44),type="l",col="red")
legend("topleft",lty=1,col=c("black","red"),bty="n",legend = c("PPRT","PrimPol"))
abline(v=22,lty=2)
text(22,60,pos=4,"8-oxo-G")
dev.off()

pdf("../Figs/error_rate_boxplot.pdf")
par(mfrow=c(1,2))
boxplot(eR[6:26,],ylab="Error rate (%)",names = NA,main="cDNA 3' end error rate")
axis(1,at=1:6,labels=c("PPRT","Mut","PPRT|8-oxo-G",
  "Mut|8-oxo-G","HIV RT","HIV PPRT"),las=2)

boxplot(eR[32:45,],ylab="Error rate (%)",names = NA,main="cDNA 5' end error rate")
axis(1,at=1:6,labels=c("PPRT","Mut","PPRT|8-oxo-G",
                       "Mut|8-oxo-G","HIV RT","HIV PPRT"),las=2)
dev.off()

### compare err rate
pv<-c(t.test(eR$HX05[6:26],eR$HX06[6:26])$p.value,
t.test(eR$HX07[6:26],eR$HX08[6:26])$p.value,
t.test(eR$HX09[6:26],eR$HX10[6:26])$p.value,
t.test(eR$HX05[32:45],eR$HX06[32:45])$p.value,
t.test(eR$HX07[32:45],eR$HX08[32:45])$p.value,
t.test(eR$HX09[32:45],eR$HX10[32:45])$p.value)
