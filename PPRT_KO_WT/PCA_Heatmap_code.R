library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library("EnhancedVolcano")
rm(list=ls())


##Get all counts

setwd("/Users/hengyi/Documents/Projects/20221109_Guzhen_PRT/Counts")

Geno<- expand_grid(Genotype=c("WT", "PRT_KD"), Growth_time=c("2h", "5h", "8h"), rep=1:3)
Species<- as.vector(unite(Geno, col="Sample", c("Genotype", "Growth_time", "rep"), sep="_"))# Need to manually input your level everytime
Species<- Species$Sample

#bk<- unite(Geno, col="bkg", c("Genotype", "Growth_time"), sep="_")
#Geno1<- bk %>% dplyr::select(-rep) %>% cbind(Geno, Species)

Hd_input<- c("Scaffold", "Start", "End", "Gene_name", "Score","Strand", "Type", "Counts")

#to get the information of the dataset number and names
files <- list.files(path= ".", pattern=".counts", ignore.case=T)

#merge
mydata<-read.table(files[1], col.names = Hd_input)
for(i in 2:length(files))
{
  data <- read.table(files[i], col.names = Hd_input)
  
  mydata <- merge(mydata, data, by= c("Scaffold", "Start", "End", "Gene_name", "Score","Strand", "Type"), all=TRUE)
  
  rm(data)
}

new_Col_names<-  c(c("Scaffold", "Start", "End", "Gene_name", "Score","Strand", "Type"), Species)
dat_all<-mydata
rm(mydata)
colnames(dat_all)<-new_Col_names
dat_all[is.na(dat_all)]<-0

SC366_genes_info<- read.csv("../SC366_gene_full_info",sep="\t")
colnames(SC366_genes_info)<- c("Scaffold", "Start", "End", "Ref_seq", "Symbol", "Locus_tag_name", "Description")

dat_all$id<- paste(dat_all$Scaffold, dat_all$Start, dat_all$End, sep="_")
SC366_genes_info$id<- paste(SC366_genes_info$Scaffold, SC366_genes_info$Start, SC366_genes_info$End, sep="_")

dat_final<- dat_all %>% full_join(SC366_genes_info, by="id") %>% select(c(1:4, 30,32,33,6:25))

setwd("../Results/Profiles/")
gene_info<- dat_final[,1:9]
colnames(gene_info)<- c("Scaffold", "Start", "End", "Gene_name", "Ref_seq", "Locus_tag_name", "Description", "Strand", "Type")

dat_temp<- dat_final[,10:27] +1
dat_1<- apply(dat_temp, 2, function(x){x/sum(x) * 1000000})
dat_CPM<-cbind(gene_info,dat_1)

#dat_CPM[dat_CPM$Type=="ncRNA", ]$Type<- "Other_ncRNA" 
#dat_CPM[dat_CPM$Type=="protein_coding", ]$Type<- "Protein_coding" 
dat_CPM$Type <- sub("ncRNA", "Other_ncRNA", dat_CPM$Type)
dat_CPM$Type <- sub("protein_coding", "Protein_coding", dat_CPM$Type)
dat_CPM$Type<- factor(dat_CPM$Type, levels = c("tRNA","tmRNA", "Protein_coding", "Other_ncRNA", "rRNA"))


#total RNA profiling
pdf("./total_Profile.pdf")
agg<-aggregate(.~Type,data=dat_CPM[dat_CPM$Type!="ERCC"&dat_CPM$Type!="SP",c(9,10:27)],FUN=sum)
agg.prop<-data.frame(prop.table(as.matrix(agg[,2:(1+length(files))]),margin = 2))
rownames(agg.prop)<-agg$Type
colnames(agg.prop)<-Species
scol <- c(brewer.pal(12, "Set3"))[c(1,2,3,4,6,5,7,8,9,10,11,12)]

mp<-barplot(as.matrix(agg.prop[, c(1:length(files))]), col=scol, axes=F, width=0.04,space=0.2,cex.names=0.6,
            ylim=c(0,1),legend.text = agg$Type,
            las=2,
            args.legend=list(x=1.25,y=0.5,bty="n",cex=1.0,yjust=0.5),
            adj=0.12,xlim=c(0,1.2)
)

axis(1,labels=c(NA,NA),las=1,at=mp,pos=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,0.25,0.5,0.75,1),pos=0)
mtext(side=3,at=mean(mp),line=0.5,text="PRT_knockout",cex=1)

dev.off()

#total RNA profiling without rRNA
pdf("./total_Profile(without_rRNA).pdf")
agg<-aggregate(.~Type,data=dat_CPM[dat_CPM$Type!="ERCC"&dat_CPM$Type!="SP"&dat_CPM$Type!="rRNA",c(9,10:27)],FUN=sum)
agg.prop<-data.frame(prop.table(as.matrix(agg[,2:(1+length(files))]),margin = 2))
rownames(agg.prop)<-agg$Type
colnames(agg.prop)<-Species
scol <- c(brewer.pal(12, "Set3"))[c(1,2,3,4,6,5,7,8,9,10,11,12)]

mp<-barplot(as.matrix(agg.prop[, c(1:length(files))]), col=scol, axes=F, width=0.04,space=0.2,cex.names=0.6,
            ylim=c(0,1),legend.text = agg$Type,
            las=2,
            args.legend=list(x=1.25,y=0.5,bty="n",cex=1.0,yjust=0.5),
            adj=0.12,xlim=c(0,1.2)
)

axis(1,labels=c(NA,NA),las=1,at=mp,pos=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,0.25,0.5,0.75,1),pos=0)
mtext(side=3,at=mean(mp),line=0.5,text="PRT_knockout (no rRNA)",cex=1)

dev.off()


##DESEQ

#Load library
#BiocManager::install("DESeq2")
library(DESeq2)
library(ashr)

coldata_DE<- as.data.frame(Geno)
rownames(coldata_DE)<- Species
coldata_DE$Group<- paste0(coldata_DE$Genotype, coldata_DE$Growth_time)
coldata_DE$Group<- factor(coldata_DE$Group, level=c("WT2h", "WT5h", "WT8h", "PRT_KD2h", "PRT_KD5h", "PRT_KD8h"))

#coldata_DE$Genotype<- factor(coldata_DE$Genotype, levels=c("WT", "PRT_KD"))
#coldata_DE$Genotype<- factor(coldata_DE$Growth_time, levels=c("2h", "5h", "8h"))


dat_DE<- dat_final[dat_final$Type!="rRNA",c(6, 10:27)]
rownames(dat_DE)<- dat_DE$Locus_tag_name
dat_DE<- dat_DE[,2:19]
dat_DE<- dat_DE + 1

dds1<- DESeqDataSetFromMatrix(countData= dat_DE,
                              colData= coldata_DE,
                              design= ~ Genotype + Growth_time + Genotype:Growth_time )

dds1<- DESeq(dds1, betaPrior = F)


#rld transform
rld1<- rlog(dds1, blind=F)

pdf("./PCA_all.pdf")
pcaData <- plotPCA(rld1, intgroup=c("Genotype", "Growth_time"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, shape=Genotype, color=Growth_time)) +
  geom_jitter(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_x_continuous(limits=c(-60, 60)) + 
  scale_y_continuous(limits=c(-60,60)) +
  theme_classic()
dev.off()

resultsNames(dds1) # this is to understand the x terms in the linear model

#KO vs WT (b1), aka, difference between KO and WT at 2h
resultsNames(dds1)[2]
res1<- results(dds1, contrast=c(0,1,0,0,0,0))
res1<- lfcShrink(dds1, contrast=c(0,1,0,0,0,0), res=res1, type='ashr')

res1$Locus_tag_name<- rownames(res1)
res1_df<- as.data.frame(res1)
res1_df<- dat_CPM[, c(1:9)] %>% right_join(res1_df, by="Locus_tag_name")

#5h vs 2h (b2), aka the difference between 5h and 2h in WT
resultsNames(dds1)[3]
res2<- results(dds1, contrast=c(0,0,1,0,0,0))
res2<- lfcShrink(dds1, contrast=c(0,0,1,0,0,0), res=res2, type='ashr')

res2$Locus_tag_name<- rownames(res2)
res2_df<- as.data.frame(res2)
res2_df<- dat_CPM[, c(1:9)] %>% right_join(res2_df, by="Locus_tag_name")

#8h vs 2h (b3), aka the difference between 8h and 2h in WT
resultsNames(dds1)[4]
res3<- results(dds1, contrast=c(0,0,0,1,0,0))
res3<- lfcShrink(dds1, contrast=c(0,0,0,1,0,0), res=res3, type='ashr')

res3$Locus_tag_name<- rownames(res3)
res3_df<- as.data.frame(res3)
res3_df<- dat_CPM[, c(1:9)] %>% right_join(res3_df, by="Locus_tag_name")

#5h_KO interaction term (b4), aka the specific change in KO at 5h 
resultsNames(dds1)[5]
res4<- results(dds1, contrast=c(0,0,0,0,1,0))
res4<- lfcShrink(dds1, contrast=c(0,0,0,0,1,0), res=res4, type='ashr')

res4$Locus_tag_name<- rownames(res4)
res4_df<- as.data.frame(res4)
res4_df<- dat_CPM[, c(1:9)] %>% right_join(res4_df, by="Locus_tag_name")

#8h_KO interaction term (b5), aka the specific change in KO at 8h 
resultsNames(dds1)[6]
res5<- results(dds1, contrast=c(0,0,0,0,0,1))
res5<- lfcShrink(dds1, contrast=c(0,0,0,0,0,1), res=res5, type='ashr')

res5$Locus_tag_name<- rownames(res5)
res5_df<- as.data.frame(res5)
res5_df<- dat_CPM[, c(1:9)] %>% right_join(res5_df, by="Locus_tag_name")

res1_sig<- res1_df %>% filter(padj< 0.001) %>% filter(abs(log2FoldChange)>2)
res2_sig<- res2_df %>% filter(padj< 0.001) %>% filter(abs(log2FoldChange)>2)
res3_sig<- res3_df %>% filter(padj< 0.001) %>% filter(abs(log2FoldChange)>2)
res4_sig<- res4_df %>% filter(padj< 0.001) %>% filter(abs(log2FoldChange)>2)
res5_sig<- res5_df %>% filter(padj< 0.001) %>% filter(abs(log2FoldChange)>2)

res1_sig_less<- res1_df %>% filter(padj< 0.05) 
res2_sig_less<- res2_df %>% filter(padj< 0.05) 
res3_sig_less<- res3_df %>% filter(padj< 0.05) 
res4_sig_less<- res4_df %>% filter(padj< 0.05) 
res5_sig_less<- res5_df %>% filter(padj< 0.05) 

gldata1<- assay(rld1)
Z1<- t(scale(t(gldata1)))
sampleDists<- dist(t(gldata1))
sampleDistMatrix<- as.matrix(sampleDists)
colnames(sampleDistMatrix)<- NULL
colors<- colorRampPalette(rev(brewer.pal(9, "Blues"))) (255)

pheatmap(sampleDistMatrix,
         cluster_rows=F,
         cluster_cols = F,
         col=colors)



#rm(list=ls())

#This needs to be changed based on login computer 
Box_path<- "/Users/hengyi/Documents/"

Folder_path<- paste0(Box_path, "Projects/20221109_Guzhen_pRT/20240620_Revisit_with_higherDEseq/ShinyGO_output")

setwd(Folder_path)

#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("plyr")
#install.packages("ggplot2")
#install.packages("RColorBrewer")
#install.packages("ggpubr")
#install.packages("rstatix")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("EnhancedVolcano")
#install.packages("pheatmap", force=T)

#install.packages("matrixStats")

library(matrixStats)
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library("EnhancedVolcano")
library(pheatmap)

#KEGG
files = list.files(path= ".", pattern="KEGG", ignore.case=T)
#Hd_input<- c("Enrichment_FDR", "nGenes", "nPathway_Genes", "Fold_Enrichment", "Pathway", "URL", "Genes")

for(i in 1:length(files)){ 
  if (i==1) {
    mydata<-read.table(files[i], sep=",", header=T, colClasses = c(NA, "NULL", "NULL", "NULL", NA, "NULL", "NULL"))
  } else {
    data <- read.table(files[i], sep=",", header=T, colClasses = c(NA, "NULL", "NULL", "NULL", NA, "NULL", "NULL"))
    
    mydata <- merge(mydata, data, by= "Pathway", all=T)
    
    rm(data)
  }
}

colnames(mydata)<- c("Pathway", "KO_vs_WT(2h)", "KO_vs_WT(5h)", "KO_vs_WT(8h)")

mydata[is.na(mydata)] <- 1

#keep<- rowMins(as.matrix(mydata[,2:6])) < 0.01

KEGG<- mydata

#BP
files = list.files(path= ".", pattern="BP", ignore.case=T)
#Hd_input<- c("Enrichment_FDR", "nGenes", "nPathway_Genes", "Fold_Enrichment", "Pathway", "URL", "Genes")

for(i in 1:length(files)){ 
  if (i==1) {
    mydata<-read.table(files[i], sep=",", header=T, colClasses = c(NA, "NULL", "NULL", "NULL", NA, "NULL", "NULL"))
  } else {
    data <- read.table(files[i], sep=",", header=T, colClasses = c(NA, "NULL", "NULL", "NULL", NA, "NULL", "NULL"))
    
    mydata <- merge(mydata, data, by= "Pathway", all=T)
    
    rm(data)
  }
}

colnames(mydata)<- c("Pathway", "KO_vs_WT(2h)", "KO_vs_WT(5h)", "KO_vs_WT(8h)")

mydata[is.na(mydata)] <- 1

#keep<- rowMins(as.matrix(mydata[,2:6])) < 0.01

BP<- mydata


#CC
files = list.files(path= ".", pattern="CC", ignore.case=T)
#Hd_input<- c("Enrichment_FDR", "nGenes", "nPathway_Genes", "Fold_Enrichment", "Pathway", "URL", "Genes")

for(i in 1:length(files)){ 
  if (i==1) {
    mydata<-read.table(files[i], sep=",", header=T, colClasses = c(NA, "NULL", "NULL", "NULL", NA, "NULL", "NULL"))
  } else {
    data <- read.table(files[i], sep=",", header=T, colClasses = c(NA, "NULL", "NULL", "NULL", NA, "NULL", "NULL"))
    
    mydata <- merge(mydata, data, by= "Pathway", all=T)
    
    rm(data)
  }
}

colnames(mydata)<- c("Pathway", "KO_vs_WT(2h)", "KO_vs_WT(5h)", "KO_vs_WT(8h)")

mydata[is.na(mydata)] <- 1

#keep<- rowMins(as.matrix(mydata[,2:6])) < 0.01

CC<- mydata

#MF
files = list.files(path= ".", pattern="MF", ignore.case=T)
#Hd_input<- c("Enrichment_FDR", "nGenes", "nPathway_Genes", "Fold_Enrichment", "Pathway", "URL", "Genes")

for(i in 1:length(files)){ 
  if (i==1) {
    mydata<-read.table(files[i], sep=",", header=T, colClasses = c(NA, "NULL", "NULL", "NULL", NA, "NULL", "NULL"))
  } else {
    data <- read.table(files[i], sep=",", header=T, colClasses = c(NA, "NULL", "NULL", "NULL", NA, "NULL", "NULL"))
    
    mydata <- merge(mydata, data, by= "Pathway", all=T)
    
    rm(data)
  }
}

colnames(mydata)<- c("Pathway", "KO_vs_WT(2h)", "KO_vs_WT(5h)", "KO_vs_WT(8h)")

mydata[is.na(mydata)] <- 1

#keep<- rowMins(as.matrix(mydata[,2:6])) < 0.01

MF<- mydata

#Ecocyc
files = list.files(path= ".", pattern="Ecocyc", ignore.case=T)
#Hd_input<- c("Enrichment_FDR", "nGenes", "nPathway_Genes", "Fold_Enrichment", "Pathway", "URL", "Genes")

for(i in 1:length(files)){ 
  if (i==1) {
    mydata<-read.table(files[i], sep=",", header=T, colClasses = c(NA, "NULL", "NULL", "NULL", NA, "NULL", "NULL"))
  } else {
    data <- read.table(files[i], sep=",", header=T, colClasses = c(NA, "NULL", "NULL", "NULL", NA, "NULL", "NULL"))
    
    mydata <- merge(mydata, data, by= "Pathway", all=T)
    
    rm(data)
  }
}

colnames(mydata)<- c("Pathway", "KO_vs_WT(2h)", "KO_vs_WT(5h)", "KO_vs_WT(8h)")

mydata[is.na(mydata)] <- 1

#keep<- rowMins(as.matrix(mydata[,2:6])) < 0.01

Ecocyc<- mydata


KEGG$GSEA<- "KEGG"
BP$GSEA<- "BP"
CC$GSEA<- "CC"
MF$GSEA<- "MF"
Ecocyc$GSEA<- "Ecocyc"

All_GSEA<- rbind(KEGG, Ecocyc, BP, CC, MF)


###################################################################
#Visualization
#########################################################
setwd("../Results/")

#Combined GSEA for all 3 cases
rownames(All_GSEA)<- All_GSEA$Pathway


#Ecocyc
df1<- All_GSEA %>% filter(GSEA=="Ecocyc") %>% select(c(2:4))
df1[is.na(df1)]<-1

df_log<- df1 %>% log10()
df_log<- 0-df_log

pdf("./Ecocyc_GSEA.pdf")
p<- pheatmap(t(as.matrix(df_log)), cluster_cols = T, cluster_rows = F, fontsize=6)
print(p)
dev.off()


#KEGG
df1<- All_GSEA %>% filter(GSEA=="KEGG") %>% select(c(2:4))
df1[is.na(df1)]<-1

df_log<- df1 %>% log10()
df_log<- 0-df_log
df_log[df_log>3]<- 3

pdf("./KEGG_GSEA.pdf")
p<- pheatmap(t(as.matrix(df_log)), cluster_cols = T, cluster_rows = F, fontsize=6)
print(p)
dev.off()


#BP
df1<- All_GSEA %>% filter(GSEA=="BP") %>% select(c(2:4))
df1[is.na(df1)]<-1

df_log<- df1 %>% log10()
df_log<- 0-df_log
df_log[df_log>4]<- 4

pdf("./BP_GSEA.pdf")
p<- pheatmap(t(as.matrix(df_log)), cluster_cols = T, cluster_rows = F, fontsize=6)
print(p)
dev.off()

#CC
df1<- All_GSEA %>% filter(GSEA=="CC") %>% select(c(2:4))
df1[is.na(df1)]<-1

df_log<- df1 %>% log10()
df_log<- 0-df_log

pdf("./CC_GSEA.pdf")
p<- pheatmap(t(as.matrix(df_log)), cluster_cols = T, cluster_rows = F, fontsize=9)
print(p)
dev.off()

#MF
df1<- All_GSEA %>% filter(GSEA=="MF") %>% select(c(2:4))
df1[is.na(df1)]<-1

df_log<- df1 %>% log10()
df_log<- 0-df_log
df_log[df_log>3]<- 3

pdf("./MF_GSEA.pdf")
p<- pheatmap(t(as.matrix(df_log)), cluster_cols = T, cluster_rows = F, fontsize=6)
print(p)
dev.off()

#Z-matrix preparation
gldata1<- assay(rld1)
Z1<- t(scale(t(gldata1)))

##5h BP genes
#install.packages("tidyverse")
library(tidyverse)
setwd("../ShinyGO_output/")

files = list.files(path= ".", pattern="BP", ignore.case=T)
p1<- read_csv(files[1])
p1_genes<- p1[,7] %>% as.vector %>% unlist()
p1_genes <- strsplit(p1_genes, " ") %>% unlist()

p2<- read_csv(files[2])
p2_genes<- p2[,7] %>% as.vector %>% unlist()
p2_genes <- strsplit(p2_genes, " ") %>% unlist()

p3<- read_csv(files[3])
p3_genes<- p3[,7] %>% as.vector %>% unlist()
p3_genes <- strsplit(p3_genes, " ") %>% unlist()

BP_genes_5h<- p2_genes %>% unique()
BP_genes_8h<- p3_genes %>% unique()

BP_5h_Z<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% BP_genes_5h)

Z1_5h_BP<- Z1[rownames(Z1) %in% BP_5h_Z$Locus_tag_name,]
Z1_5h_BP<- as.data.frame(Z1_5h_BP)
Z1_5h_BP$Locus_tag_name<- rownames(Z1_5h_BP)
Z1_5h_BP_df<- Z1_5h_BP %>% left_join(BP_5h_Z, by="Locus_tag_name")
Z1_5h_BP_df$label<- paste(Z1_5h_BP_df$Gene_name, Z1_5h_BP_df$Locus_tag_name, Z1_5h_BP_df$Description, sep = ":")
rownames(Z1_5h_BP_df)<- Z1_5h_BP_df$label
Z1_5h_BP_m<- Z1_5h_BP_df[,1:18] %>% as.matrix() 
p<- pheatmap(Z1_5h_BP_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=6)

###Selection of 5h specific changes: the DNA damage, oxidative stress, anaerobic respiration, etc
##BP
files = list.files(path= ".", pattern="BP", ignore.case=T)
p<- read_csv(files[2])

p_genes<- p %>% filter(Pathway=="GO:0006974 cellular response to DNA damage stimulus") %>% select(7) %>% as.vector %>% unlist()
p_genes <- strsplit(p_genes, " ") %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% p_genes)
Z1_p<- Z1[rownames(Z1) %in% gene_info$Locus_tag_name,]
Z1_p<- as.data.frame(Z1_p)
Z1_p$Locus_tag_name<- rownames(Z1_p)
Z1_p_df<- Z1_p %>% left_join(gene_info, by="Locus_tag_name")
Z1_p_df$label<- paste(Z1_p_df$Gene_name, Z1_p_df$Locus_tag_name, Z1_p_df$Description, sep = ":")
rownames(Z1_p_df)<- Z1_p_df$label
Z1_p_m<- Z1_p_df[,1:18] %>% as.matrix()
Z1_p_m[Z1_p_m < -2] <- -2
Z1_p_m[Z1_p_m > 2] <- 2
q<- pheatmap(Z1_p_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=6)
pdf("../Results/GO_0006974 cellular response to DNA damage stimulus.pdf")
print(q)
dev.off()

p_genes<- p %>% filter(Pathway=="GO:0009061 anaerobic respiration") %>% select(7) %>% as.vector %>% unlist()
p_genes <- strsplit(p_genes, " ") %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% p_genes)
Z1_p<- Z1[rownames(Z1) %in% gene_info$Locus_tag_name,]
Z1_p<- as.data.frame(Z1_p)
Z1_p$Locus_tag_name<- rownames(Z1_p)
Z1_p_df<- Z1_p %>% left_join(gene_info, by="Locus_tag_name")
Z1_p_df$label<- paste(Z1_p_df$Gene_name, Z1_p_df$Locus_tag_name, Z1_p_df$Description, sep = ":")
rownames(Z1_p_df)<- Z1_p_df$label
Z1_p_m<- Z1_p_df[,1:18] %>% as.matrix()
Z1_p_m[Z1_p_m < -2] <- -2
Z1_p_m[Z1_p_m > 2] <- 2
q<- pheatmap(Z1_p_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=6)
pdf("../Results/GO_0009061 anaerobic respiration.pdf")
print(q)
dev.off()

p_genes<- p %>% filter(Pathway=="GO:0019646 aerobic electron transport chain") %>% select(7) %>% as.vector %>% unlist()
p_genes <- strsplit(p_genes, " ") %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% p_genes)
Z1_p<- Z1[rownames(Z1) %in% gene_info$Locus_tag_name,]
Z1_p<- as.data.frame(Z1_p)
Z1_p$Locus_tag_name<- rownames(Z1_p)
Z1_p_df<- Z1_p %>% left_join(gene_info, by="Locus_tag_name")
Z1_p_df$label<- paste(Z1_p_df$Gene_name, Z1_p_df$Locus_tag_name, Z1_p_df$Description, sep = ":")
rownames(Z1_p_df)<- Z1_p_df$label
Z1_p_m<- Z1_p_df[,1:18] %>% as.matrix() 
Z1_p_m[Z1_p_m < -2] <- -2
Z1_p_m[Z1_p_m > 2] <- 2
q<- pheatmap(Z1_p_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=6)
pdf("../Results/GO_0019646 aerobic electron transport chain.pdf")
print(q)
dev.off()

p_genes<- p %>% filter(Pathway=="GO:0070814 hydrogen sulfide biosynthetic process") %>% select(7) %>% as.vector %>% unlist()
p_genes <- strsplit(p_genes, " ") %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% p_genes)
Z1_p<- Z1[rownames(Z1) %in% gene_info$Locus_tag_name,]
Z1_p<- as.data.frame(Z1_p)
Z1_p$Locus_tag_name<- rownames(Z1_p)
Z1_p_df<- Z1_p %>% left_join(gene_info, by="Locus_tag_name")
Z1_p_df$label<- paste(Z1_p_df$Gene_name, Z1_p_df$Locus_tag_name, Z1_p_df$Description, sep = ":")
rownames(Z1_p_df)<- Z1_p_df$label
Z1_p_m<- Z1_p_df[,1:18] %>% as.matrix()
Z1_p_m[Z1_p_m < -2] <- -2
Z1_p_m[Z1_p_m > 2] <- 2
q<- pheatmap(Z1_p_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=6)
pdf("../Results/GO_0070814 hydrogen sulfide biosynthetic process.pdf")
print(q)
dev.off()

p_genes<- p %>% filter(Pathway=="GO:0008272 sulfate transport") %>% select(7) %>% as.vector %>% unlist()
p_genes <- strsplit(p_genes, " ") %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% p_genes)
Z1_p<- Z1[rownames(Z1) %in% gene_info$Locus_tag_name,]
Z1_p<- as.data.frame(Z1_p)
Z1_p$Locus_tag_name<- rownames(Z1_p)
Z1_p_df<- Z1_p %>% left_join(gene_info, by="Locus_tag_name")
Z1_p_df$label<- paste(Z1_p_df$Gene_name, Z1_p_df$Locus_tag_name, Z1_p_df$Description, sep = ":")
rownames(Z1_p_df)<- Z1_p_df$label
Z1_p_m<- Z1_p_df[,1:18] %>% as.matrix()
Z1_p_m[Z1_p_m < -2] <- -2
Z1_p_m[Z1_p_m > 2] <- 2
q<- pheatmap(Z1_p_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=6)
pdf("../Results/GO_0008272 sulfate transport.pdf")
print(q)
dev.off()

p_genes<- p %>% filter(Pathway=="GO:0006979 response to oxidative stress") %>% select(7) %>% as.vector %>% unlist()
p_genes <- strsplit(p_genes, " ") %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% p_genes)
Z1_p<- Z1[rownames(Z1) %in% gene_info$Locus_tag_name,]
Z1_p<- as.data.frame(Z1_p)
Z1_p$Locus_tag_name<- rownames(Z1_p)
Z1_p_df<- Z1_p %>% left_join(gene_info, by="Locus_tag_name")
Z1_p_df$label<- paste(Z1_p_df$Gene_name, Z1_p_df$Locus_tag_name, Z1_p_df$Description, sep = ":")
rownames(Z1_p_df)<- Z1_p_df$label
Z1_p_m<- Z1_p_df[,1:18] %>% as.matrix()
Z1_p_m[Z1_p_m < -2] <- -2
Z1_p_m[Z1_p_m > 2] <- 2
q<- pheatmap(Z1_p_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=6)
pdf("../Results/GO_0006979 response to oxidative stress.pdf")
print(q)
dev.off()

##Ecocyc
files = list.files(path= ".", pattern="Ecocyc", ignore.case=T)
p<- read_csv(files[2])

p_genes<- p %>% filter(Pathway=="PWY0-1576 hydrogen to fumarate electron transfer") %>% select(7) %>% as.vector %>% unlist()
p_genes <- strsplit(p_genes, " ") %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% p_genes)
Z1_p<- Z1[rownames(Z1) %in% gene_info$Locus_tag_name,]
Z1_p<- as.data.frame(Z1_p)
Z1_p$Locus_tag_name<- rownames(Z1_p)
Z1_p_df<- Z1_p %>% left_join(gene_info, by="Locus_tag_name")
Z1_p_df$label<- paste(Z1_p_df$Gene_name, Z1_p_df$Locus_tag_name, Z1_p_df$Description, sep = ":")
rownames(Z1_p_df)<- Z1_p_df$label
Z1_p_m<- Z1_p_df[,1:18] %>% as.matrix() 
q<- pheatmap(Z1_p_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=6)
pdf("../Results/PWY0-1576 hydrogen to fumarate electron transfer.pdf")
print(q)
dev.off()

p_genes<- p %>% filter(Pathway=="PWY0-1353 succinate to cytochrome bd oxidase electron transfer") %>% select(7) %>% as.vector %>% unlist()
p_genes <- strsplit(p_genes, " ") %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% p_genes)
Z1_p<- Z1[rownames(Z1) %in% gene_info$Locus_tag_name,]
Z1_p<- as.data.frame(Z1_p)
Z1_p$Locus_tag_name<- rownames(Z1_p)
Z1_p_df<- Z1_p %>% left_join(gene_info, by="Locus_tag_name")
Z1_p_df$label<- paste(Z1_p_df$Gene_name, Z1_p_df$Locus_tag_name, Z1_p_df$Description, sep = ":")
rownames(Z1_p_df)<- Z1_p_df$label
Z1_p_m<- Z1_p_df[,1:18] %>% as.matrix() 
q<- pheatmap(Z1_p_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=8)
pdf("../Results/PWY0-1353 succinate to cytochrome bd oxidase electron transfer.pdf")
print(q)
dev.off()

p_genes<- p %>% filter(Pathway=="PWY0-1585 formate to nitrite electron transfer") %>% select(7) %>% as.vector %>% unlist()
p_genes <- strsplit(p_genes, " ") %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% p_genes)
Z1_p<- Z1[rownames(Z1) %in% gene_info$Locus_tag_name,]
Z1_p<- as.data.frame(Z1_p)
Z1_p$Locus_tag_name<- rownames(Z1_p)
Z1_p_df<- Z1_p %>% left_join(gene_info, by="Locus_tag_name")
Z1_p_df$label<- paste(Z1_p_df$Gene_name, Z1_p_df$Locus_tag_name, Z1_p_df$Description, sep = ":")
rownames(Z1_p_df)<- Z1_p_df$label
Z1_p_m<- Z1_p_df[,1:18] %>% as.matrix() 
q<- pheatmap(Z1_p_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=6)
pdf("../Results/PWY0-1585 formate to nitrite electron transfer.pdf")
print(q)
dev.off()

##MF
files = list.files(path= ".", pattern="MF", ignore.case=T)
p<- read_csv(files[2])

p_genes<- p %>% filter(Pathway=="GO:0016491 oxidoreductase activity") %>% select(7) %>% as.vector %>% unlist()
p_genes <- strsplit(p_genes, " ") %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Gene_name %in% p_genes)
Z1_p<- Z1[rownames(Z1) %in% gene_info$Locus_tag_name,]
Z1_p<- as.data.frame(Z1_p)
Z1_p$Locus_tag_name<- rownames(Z1_p)
Z1_p_df<- Z1_p %>% left_join(gene_info, by="Locus_tag_name")
Z1_p_df$label<- paste(Z1_p_df$Gene_name, Z1_p_df$Locus_tag_name, Z1_p_df$Description, sep = ":")
rownames(Z1_p_df)<- Z1_p_df$label
Z1_p_m<- Z1_p_df[,1:18] %>% as.matrix() 
Z1_p_m[Z1_p_m < -2] <- -2
Z1_p_m[Z1_p_m > 2] <- 2
q<- pheatmap(Z1_p_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=6)
pdf("../Results/GO_0016491 oxidoreductase activity.pdf")
print(q)
dev.off()

ncRNA_genes<- dat_CPM %>% filter(Type=="Other_ncRNA") %>% select(Locus_tag_name) %>% as.vector() %>% unlist()
gene_info<- dat_final[,c(4:7)] %>% filter(dat_final$Locus_tag_name %in% ncRNA_genes)
Z1_nc<- Z1[rownames(Z1) %in% ncRNA_genes,]
Z1_nc<- as.data.frame(Z1_nc)
Z1_nc$Locus_tag_name<- rownames(Z1_nc)
Z1_nc_df<- Z1_nc %>% left_join(gene_info, by="Locus_tag_name")
#Z1_nc_df$label<- paste(Z1_nc_df$Gene_name, Z1_nc_df$Locus_tag_name, Z1_nc_df$Description, sep = ":")
rownames(Z1_nc_df)<- Z1_nc_df$Gene_name
Z1_nc_m<- Z1_nc_df[,1:18] %>% as.matrix()
Z1_nc_m[Z1_nc_m < -2] <- -2
Z1_nc_m[Z1_nc_m > 2] <- 2
q<- pheatmap(Z1_nc_m, annotation_col = coldata_DE[,1:2], cluster_cols = F, cluster_rows = T, fontsize=10)
pdf("../Results/ncRNA_heatmap.pdf")
print(q)
dev.off()

