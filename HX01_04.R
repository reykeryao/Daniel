rm(list=ls())
library(tidyverse)
library(RColorBrewer)
setwd("/stor/work/Lambowitz/yaojun/Work/JA25159_Daniel/")
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
HX01<-read.delim("Merged/HX01.info")
HX02<-read.delim("Merged/HX02.info")
HX03<-read.delim("Merged/HX03.info")
HX04<-read.delim("Merged/HX04.info")
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