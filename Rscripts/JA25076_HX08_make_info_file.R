library(stringr)
dat<-read.table(gzfile(paste0("Daniel/Merged/HX08_JA25076.tab.gz")),col.names=c("ID","Seq"))
dat$Len<-nchar(dat$Seq)
dat$A<-str_count(dat$Seq,"A")
dat$C<-str_count(dat$Seq,"C")
dat$G<-str_count(dat$Seq,"G")
dat$T<-str_count(dat$Seq,"T")
dat$polyA<-0
dat$polyC<-0
dat$polyG<-0
dat$polyT<-0
for (i in 1:length(dat$Seq)){
 if (i%%10000==0){print(i)}
 seq<-dat$Seq[i]
 seq<-unique(unlist(str_extract_all(seq, "(.)\\1+")))
 polyA<-seq[grep("A",seq)]
 polyC<-seq[grep("C",seq)]
 polyG<-seq[grep("G",seq)]
 polyT<-seq[grep("T",seq)] 
 dat$polyA[i]<-max(nchar(polyA))
 dat$polyC[i]<-max(nchar(polyC))
 dat$polyG[i]<-max(nchar(polyG))
 dat$polyT[i]<-max(nchar(polyT))
 }

gz1 <- gzfile(paste0("Daniel/Merged/HX08_JA25076.info.gz"), "w")
write.table(dat,gz1,sep="\t",row.names=F,quote=F)  
close(gz1)

