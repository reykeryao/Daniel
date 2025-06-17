rm(list=ls())
library(tidyverse)
library(RColorBrewer)
library(Biostrings)
setwd("/stor/work/Lambowitz/yaojun/Work/JA25159_Daniel/Daniel/")
### overall NT frq
pdf("Figs/Nt_frq_total.pdf")
frq<-read.delim("nt_frq.txt",row.names = 1)
frq<-t(frq)[1:4,c(2,1,4,3)]
barplot(frq,names.arg = c("WT","WT+dNTP","Mut","Mut+dNTP"),
        col=c("lightblue","blue","gold","tomato"),ylim=c(0,3e7),
        ylab="Reads",main="Nucleotide frequency",yaxt="n")
legend("topleft",legend = rownames(frq),fill=c("lightblue","blue","gold","tomato"),bty="n")
axis(2,at=c(0,1e7,2e7,3e7),labels = expression(0, 10^7,2*x*10^7,3*x*10^7),las=2)
dev.off()

### polyN distribution
polyN<-read.table("polyN.txt",row.names = 1,header = TRUE)
pdf("Figs/polyN.pdf")
par(mfrow=c(2,2))
##polyA
pick=c(2:10)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads contain polyA")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
legend("topright",bty="n",legend = c("WT","WT+dNTP","Mut","Mut+dNTP"),
       lty=c(2,1,2,1),lwd=c(1,1.5,1,1.5),col=c("black","black","red","red"))
##polyC
pick=c(11:19)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads contain polyC")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
##polyG
pick=c(20:28)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads contain polyG")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)

##polyT
pick=c(29:37)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads contain polyT")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
dev.off()

### density of NT in each reads
HX01<-read.delim(gzfile("Merged/HX01.info.gz"))
HX02<-read.delim(gzfile("Merged/HX02.info.gz"))
HX03<-read.delim(gzfile("Merged/HX03.info.gz"))
HX04<-read.delim(gzfile("Merged/HX04.info.gz"))
pdf("Figs/nt_frq_density.pdf")
par(mfrow=c(2,2))
plot(density(100*HX02$A/HX02$Len),type="l",xlim=c(0,100),main="WT",ylim=c(0,0.2),
     lwd=2,col="lightblue",xlab="Nucleotide frequency (%)")
lines(density(100*HX02$C/HX02$Len),col="blue",lwd=2)
lines(density(100*HX02$G/HX02$Len),col="gold",lwd=2)
lines(density(100*HX02$T/HX02$Len),col="tomato",lwd=2)
abline(v=25,lty=2)
legend("topright",legend = c("A","C","G","T"),col=c("lightblue","blue","gold","tomato"),
       lwd=2,bty="n")

plot(density(100*HX01$A/HX01$Len),type="l",xlim=c(0,100),main="WT+dNTP",ylim=c(0,0.1),
     lwd=2,col="lightblue",xlab="Nucleotide frequency (%)")
lines(density(100*HX01$C/HX01$Len),col="blue",lwd=2)
lines(density(100*HX01$G/HX01$Len),col="gold",lwd=2)
lines(density(100*HX01$T/HX01$Len),col="tomato",lwd=2)
abline(v=25,lty=2)

plot(density(100*HX04$A/HX04$Len),type="l",xlim=c(0,100),main="Mut",ylim=c(0,0.15),
     lwd=2,col="lightblue",xlab="Nucleotide frequency (%)")
lines(density(100*HX04$C/HX04$Len),col="blue",lwd=2)
lines(density(100*HX04$G/HX04$Len),col="gold",lwd=2)
lines(density(100*HX04$T/HX04$Len),col="tomato",lwd=2)
abline(v=25,lty=2)

plot(density(100*HX03$A/HX03$Len),type="l",xlim=c(0,100),main="Mut+dNTP",ylim=c(0,0.15),
     lwd=2,col="lightblue",xlab="Nucleotide frequency (%)")
lines(density(100*HX03$C/HX03$Len),col="blue",lwd=2)
lines(density(100*HX03$G/HX03$Len),col="gold",lwd=2)
lines(density(100*HX03$T/HX03$Len),col="tomato",lwd=2)
abline(v=25,lty=2)
dev.off()

### length of polyN distribution
HX01[HX01== -Inf]<-0
HX02[HX02== -Inf]<-0
HX03[HX03== -Inf]<-0
HX04[HX04== -Inf]<-0
pdf("Figs/longest_polyA_distribution.pdf")
par(mfrow=c(2,2))
plot(density(HX02$polyA[HX02$polyA>1]),type="l",main="PolyA",xlim=c(0,20),ylim=c(0,0.5),
     lty=2,xlab="Length (nt)")
lines(density(HX01$polyA[HX01$polyA>1]),lwd=2)
lines(density(HX04$polyA[HX04$polyA>1]),col="red",lty=2)
lines(density(HX03$polyA[HX03$polyA>1]),col="red",lwd=2)
legend("topright",legend = c("WT","WT+dNTP","Mut","Mut+dNTP"),
       lty=c(2,1,2,1),lwd=c(1,1.5,1,1.5),col=c("black","black","red","red"),bty="n")

plot(density(HX02$polyC[HX02$polyC>1]),type="l",main="PolyC",xlim=c(0,20),ylim=c(0,0.5),
     lty=2,xlab="Length (nt)")
lines(density(HX01$polyC[HX01$polyC>1]),lwd=2)
lines(density(HX04$polyC[HX04$polyC>1]),col="red",lty=2)
lines(density(HX03$polyC[HX03$polyC>1]),col="red",lwd=2)
plot(density(HX02$polyG[HX02$polyG>1]),type="l",main="PolyG",xlim=c(0,20),ylim=c(0,0.5),
     lty=2,xlab="Length (nt)")
lines(density(HX01$polyG[HX01$polyG>1]),lwd=2)
lines(density(HX04$polyG[HX04$polyG>1]),col="red",lty=2)
lines(density(HX03$polyG[HX03$polyG>1]),col="red",lwd=2)
plot(density(HX02$polyT[HX02$polyT>1]),type="l",main="PolyT",xlim=c(0,20),ylim=c(0,0.5),
     lty=2,xlab="Length (nt)")
lines(density(HX01$polyT[HX01$polyT>1]),lwd=2)
lines(density(HX04$polyT[HX04$polyT>1]),col="red",lty=2)
lines(density(HX03$polyT[HX03$polyT>1]),col="red",lwd=2)
dev.off()

### reads start with polyN
polyN<-read.table("start_polyN.txt",row.names = 1,header = TRUE)
pdf("Figs/start_polyN.pdf")
par(mfrow=c(2,2))
##polyA
pick=c(2:10)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads start with polyA")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
legend("topright",bty="n",legend = c("WT","WT+dNTP","Mut","Mut+dNTP"),
       lty=c(2,1,2,1),lwd=c(1,1.5,1,1.5),col=c("black","black","red","red"))
##polyC
pick=c(11:19)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads start with polyC")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
##polyG
pick=c(20:28)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads start with polyG")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)

##polyT
pick=c(29:37)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads start with polyT")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
dev.off()


### reads end with polyN
polyN<-read.table("end_polyN.txt",row.names = 1,header = TRUE)
pdf("Figs/end_polyN.pdf")
par(mfrow=c(2,2))
##polyA
pick=c(2:10)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads end with polyA")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
legend("topright",bty="n",legend = c("WT","WT+dNTP","Mut","Mut+dNTP"),
       lty=c(2,1,2,1),lwd=c(1,1.5,1,1.5),col=c("black","black","red","red"))
##polyC
pick=c(11:19)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads end with polyC")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
##polyG
pick=c(20:28)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads end with polyG")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)

##polyT
pick=c(29:37)
plot(y=100*polyN[2,pick]/polyN$Total[2],x=c(2:10),type="l",ylim=c(0,100),
     axes=F,xlab="Polynucleotide length (nt)",ylab="Reads (%)",lty=2,
     main="Reads end with polyT")
lines(y=100*polyN[1,pick]/polyN$Total[1],x=c(2:10),lwd=1.5)
lines(y=100*polyN[4,pick]/polyN$Total[4],x=c(2:10),lty=2,col="red")
lines(y=100*polyN[3,pick]/polyN$Total[3],x=c(2:10),lwd=1.5,col="red")
axis(1,at=c(2:10),labels = c(2:10))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
dev.off()

### reads polyN patterns
### polyN-NNNNNN(0 or more times)-polyN
polyN<-read.table("polyN_patterns.txt",row.names = 1,header = TRUE)
### di nucleotide (AC etc) patterns are reversed -v, need subtract form the total
for (i in 2:11){
  if (i<8){
    polyN[,i]<-(polyN[,1]-polyN[,i])/polyN[,1]
  } else {
    polyN[,i]<-polyN[,i]/polyN[,1]
  }
}

polyN<-polyN[,-1]*100
## reorgnize row
polyN<-polyN[c(2,1,4,3),]
## patterns
patterns<-paste(colnames(polyN),c("AAAANNNNCCCC or CCCCNNNNAAAA",
            "AAAANNNNGGGG or GGGGNNNNAAAA",
            "AAAANNNNTTTT or TTTTNNNNAAAA",
            "CCCCNNNNGGGG or GGGGNNNNCCCC",
            "CCCCNNNNTTTT or TTTTNNNNCCCC",
            "GGGGNNNNTTTT or TTTTNNNNGGGG",
            "AAAANNNNAAAA",
            "CCCCNNNNCCCC",
            "GGGGNNNNGGGG",
            "TTTTNNNNTTTT"),sep=": ")
  
pdf("Figs/polyN_patterns.pdf",width=8,height=11)
par(mfrow=c(2,1))
bcol<-c("lightblue","blue","gold","tomato")
mp<-barplot(as.matrix(polyN),beside = TRUE,col=bcol,axes=FALSE,ylim=c(0,100),ylab="Reads (%)")
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
legend("topleft",legend = c("WT","WT+dNTP","Mut","Mut+dNTP"),bty="n",
       pt.cex=1.5,pt.bg = bcol,pch=22,cex=0.5)

legend(x=mp[15],y=90,legend = patterns,bty="n",title = "PolyN patterns",title.adj = 0.5,cex=0.5)

bcol<-brewer.pal(10,"Set3")
mp<-barplot(t(as.matrix(polyN)),beside = TRUE,col=bcol,axes=FALSE,ylim=c(0,100),ylab="Reads (%)",
            names.arg = rep(NA,4))
axis(2,at=seq(0,100,25),labels = seq(0,100,25),las=2)
par(xpd=T)
text(colMeans(mp),100,labels = c("WT","WT+dNTP","Mut","Mut+dNTP"))
par(xpd=F)
legend(x=mp[1],y=95,legend = patterns,bty="n",title = "PolyN patterns",title.adj = 0.5, pt.cex=1.5,pt.bg = bcol,pch=22,cex=0.5)
dev.off()

### find palindrom sequence
sample<-paste0("HX0",1:4)
'''
for (i in 1:length(sample)){
  dat<-read.delim(gzfile(paste0("Merged/",sample[i],".tab.gz")),col.names=c("ID","Seq"))
  dat$left<-""
  for (j in 1:length(dat$Seq)){
    DNA<-DNAString(dat$Seq[j])
    Parl<-findPalindromes(DNA,min.armlength = 4,max.looplength = 6)
    if (length(Parl)>0){
      Left<-as.character(palindromeLeftArm(Parl))
      Left<-Left[which(nchar(Left)==max(nchar(Left)))][1]
      dat$left[j]<-Left
    } 
  }
gz1 <- gzfile(paste0("Merged/",sample[i],".parlindrom.info.gz"), "w")
write.table(dat,gz1,sep="\t",row.names=F,quote=F)  
close(gz1)
}
'''
pDict<-DNAStringSet(c("AAAA","CCCC","GGGG","TTTT"))
for (i in 1:length(sample)){
  dat<-read.delim(gzfile(paste0("Merged/",sample[i],".parlindrom.info.gz")))
  Arm<-round(sum(dat$left!="")/length(dat$left),4)
  Left<-DNAStringSet(dat$left[dat$left!=""])
  pCount<-vcountPDict(pDict,Left)
  res<-data.frame("polyN"=pDict, "Counts"=apply(pCount,1,function(x){sum(x>0)}))
  res<-rbind(res,c("Parlin",100*Arm))
  AT<-round(letterFrequency(unlist(Left),"AT")/letterFrequency(unlist(Left),"ATCG"),4)
  res<-rbind(res,c("AT",100*AT))
  if(i==1){
    tmp<-res
  } else {
    tmp<-cbind(tmp,res)
  }
}  
tmp<-tmp[,c(1,2,4,6,8)]
colnames(tmp[-1])<-sample
write.table(tmp,"Palindrome.results.txt",sep="\t",row.names=F,quote=F)
  
  
  