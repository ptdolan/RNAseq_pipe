########################
########################
#   Sleuth Analysis 
########################
########################

installPKGS<-function(){
  source("https://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
  biocLite("devtools")
  biocLite("rhdf5")
  biocLite("limma")
  devtools::install_github("stephenturner/annotables")
  devtools::install_github("pachterlab/sleuth")
  biocLite("Biostrings")
  install.packages(c("dplyr","reshape2","ggplot2","data.table","biomartr"))
  
}

annotation<-annotables::grcm38_tx2gene
DF<-data.frame(annotation)
colnames(DF)<-c("target_id","ensgene")

library(dplyr)
library(ggstance)
library(sleuth)
library(reshape2)
library(ggplot2)
library(ggpubr)

library(data.table)
library(biomaRt)
library(biomartr)
library(limma)

#######################################################################
#   BioMart - (un)comment depending on annotation type
#######################################################################

# mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",verbose = T,
#                          dataset = "mmusculus_gene_ensembl",
#                          host = "uswest.ensembl.org")
# 
# bmIDs <- biomaRt::getBM(
#   attributes = c("ensembl_gene_name","external_gene_name","refseq_mrna"),
#   mart = mart)

#######################################################################
#   SLEUTH
#######################################################################

kalDirRoot<-"/Users/ptdolan/Desktop/ENSout/"

kalDirs=list.dirs(kalDirRoot,full.names = T)

age=strsplit2(kalDirs,"_")[-c(1,14),3]
exp=strsplit2(kalDirs,"_")[-c(1,14),4]
diff<-strsplit2(exp,"-")[,2]

cond<-paste(age,diff)

input<-data.table(sample=paste(age,exp,sep = "-"),
                  age=age,
                  diff=diff,
                  path=kalDirs[c(-1,-14)])

# Initialize sleuth object
SO<- sleuth_prep(input,target_mapping = DF,gene_mode=T,
                 transformation_function = function(x) log2(x + 0.5),
                 aggregation_column = "ensgene",
                 read_bootstrap_tpm=T,extra_bootstrap_summary=T)

# Fit models
SO <- sleuth_fit(SO, ~age+diff, 'full' )
SO <- sleuth_fit(SO, ~age     , 'age'  )
SO <- sleuth_fit(SO, ~diff    , 'diff' )
SO <- sleuth_fit(SO, ~1       , 'null' )

# Compute Likelihood Ratios
SO<- sleuth_lrt( SO, null_model =  'age',  alt_model =  'full')
SO<- sleuth_lrt( SO, null_model =  'diff', alt_model =  'full')
SO<- sleuth_lrt( SO, null_model =  'null', alt_model =  'full')

SO<- sleuth_wt(SO,'age3mo')
SO<- sleuth_wt(SO,'diffD')

sleuth::sleuth_save(SO,file = "~/GitHub/RNAseq_pipe/mNSC_GeneLevel.sleuth") #"Deploys" Sleuth object as app for live viewing with shiny
sleuth::sleuth_deploy(SO,base_dir = "~/GitHub/RNAseq_pipe/")

#LR test
Diff_sleuth_table <- sleuth_results(SO, 'age:full', 'lrt', show_all = FALSE)
Age_sleuth_table <- sleuth_results(SO, 'diff:full', 'lrt', show_all = FALSE)
AgeDiff_sleuth_table <- sleuth_results(SO, 'null:full', 'lrt', show_all = FALSE)

Diff_sleuth_table<-merge(Diff_sleuth_table,annotables::grcm38,by.x = "target_id",by.y="ensgene")
Age_sleuth_table<-merge(Age_sleuth_table,  annotables::grcm38,by.x = "target_id",by.y="ensgene")

Diff_sleuth_significant <- dplyr::filter(Diff_sleuth_table, qval <= 0.01)
Age_sleuth_significant <- dplyr::filter(Age_sleuth_table, qval <= 0.01)

write.csv(Diff_sleuth_significant,file = "mNSC_Diff_Sig-q01.csv")
write.csv(Age_sleuth_significant,file = "mNSC_AgeDiff_Sig-q01.csv")

#Wald test
Wald_Diff_sleuth_table <- sleuth_results(SO,test = "diffD",test_type ='wt', show_all = F)
Wald_Age_sleuth_table <- sleuth_results(SO,test = "age3mo",test_type ='wt', show_all = F)

WDiff_sleuth_table<-merge(Wald_Diff_sleuth_table,annotables::grcm38,by.x = "target_id",by.y="ensgene")
WAge_sleuth_table<-merge(Wald_Age_sleuth_table  ,annotables::grcm38,by.x = "target_id",by.y="ensgene")

WDiff_sleuth_significant<-dplyr::filter(WDiff_sleuth_table, qval <= 0.01)
WAge_sleuth_significant<-dplyr::filter(WAge_sleuth_table, qval <= 0.01)

write.csv(WDiff_sleuth_significant,file = "mNSC_Diff_Sig-Wald-q01.csv")
write.csv(WAge_sleuth_significant,file = "mNSC_Age_Sig-Wald-q01.csv")

ggplot(SO$tests$wt$full$diffD)+geom_point(aes(b,qval,col=qval<0.01))+scale_y_log10()
ggplot(SO$tests$wt$full$age3mo)+geom_point(aes(b,qval,col=qval<0.01))+scale_y_log10()

ggplot(SO$tests$lrt$`age:full`)+geom_point(aes(pval,qval,col=qval<0.01))+scale_y_log10()
ggplot(SO$tests$lrt$`diff:full`)+geom_point(aes(pval,qval,col=qval<0.01))+scale_y_log10()

########################################################################
# PHN analysis 
########################################################################
Chaps<-read.delim("~/Downloads/Chaperones.txt")
Chaps$Names<-strsplit2(Chaps$Gene.names,split = " ")
Chaps<-reshape2::melt(Chaps)
Proteasome<-read.delim("~/Downloads/proteasome.txt")
PSMlist<-data.frame(stringsAsFactors = F,symbol=unique(c(strsplit2(Proteasome$Gene.names,split = " "))),subset="",set="Proteasome")
ChapList<-unique(data.frame(stringsAsFactors = F,symbol=Chaps$Names,subset=Chaps$Class,set="Chaperone"))

UBI<-read.delim("~/Downloads/Ubl.txt")
UBIlist<-data.frame(stringsAsFactors = F,symbol=unique(c(strsplit2(UBI$Gene.names,split = " "))),subset="",set="Ubiquitin")

allNames<-rbind.data.frame(ChapList, PSMlist, UBIlist)

modelAnnot<-function(model){
  model$ensgene<-model$target_id
  genesDF<-merge.data.frame(model,data.frame(stringsAsFactors = F,annotables::grcm38_tx2gene),by="ensgene")
  annotDF<-merge(genesDF,annotables::grcm38,by="ensgene")
  annot.DF<-merge(annotDF,allNames,by="symbol",all.x = T)
  annot.DF$qval[annot.DF$qval<1E-25]<-1e-25
  return(annot.DF)
  print("values limited to 1E-25")
}

model=SO$tests$lrt$`age:full`
LRTannot<-modelAnnot(model)

model=SO$tests$wt$full$diffD
WTannot<-modelAnnot(model)

ggplot(WTannot)+
  geom_pointrange(aes(qval,b,ymin=(b-se_b),ymax=(b+se_b),col=set))+scale_x_log10()

ggplot(WTannot[!is.na(WTannot$set),])+geom_hline(yintercept=0)+geom_vline(xintercept=25, lty=2)+
  geom_pointrange(aes(-log10(qval),b,ymin=(b-se_b),ymax=b+se_b,alpha=qval<0.01,col=subset))+
  geom_text(data=WTannot[!is.na(WTannot$set)&WTannot$qval<0.01,],aes(-log10(qval),b,label=symbol))+
  facet_wrap(~subset)

ggsave(file="~/GitHub/RNAseq_pipe/all.pdf",width=6,height=5,
ggplot(WTannot[!is.na(WTannot$set),])+
  geom_pointrangeh(aes(b,-log10(qval),xmin=(b-se_b),xmax=b+se_b,alpha=qval<0.01,col=set))+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=1e-25,lty=2)+
  geom_text(data=WTannot[!is.na(WTannot$set)&WTannot$qval<0.01,],aes(b,-log10(qval),label=symbol))
)

ggsave(file="~/GitHub/RNAseq_pipe/PSM.pdf",width=6,height=5,
ggplot(WTannot[WTannot$set=="Proteasome"&(!is.na(WTannot$set)),])+
  geom_pointrangeh(aes(b,-log10(qval),xmin=(b-se_b),xmax=b+se_b,alpha=qval<0.01,col=set))+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=1e-25,lty=2)+
  geom_text(data=WTannot[WTannot$set=="Proteasome"&WTannot$qval<0.01,],aes(b,-log10(qval),label=symbol))
)

ggsave(file="~/GitHub/RNAseq_pipe/Chaps_grid.pdf",height=8,width=4,
ggplot(WTannot[WTannot$set=="Chaperone"&(!is.na(WTannot$set))&(WTannot$subset!=""),])+
  geom_pointrangeh(aes(b,-log10(qval),xmin=(b-se_b),xmax=b+se_b,alpha=qval<0.01,col=subset))+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=1e-25,lty=2)+
  geom_text(nudge_x = .3,cex=1.8,data=unique(WTannot[WTannot$set=="Chaperone"&(!is.na(WTannot$set))&(WTannot$subset!=""),-19]),aes(b,-log10(qval),label=symbol))+
  facet_grid(subset~.,scales = 'free_y')+scale_color_brewer(palette = "Dark2")+xlim(-4,NA)
)

ggsave(file="~/GitHub/RNAseq_pipe/Chaps.pdf",height=4,width=5,
       ggplot(WTannot[WTannot$set=="Chaperone"&(!is.na(WTannot$set))&(WTannot$subset!=""),])+
         geom_pointrangeh(aes(b,-log10(qval),xmin=(b-se_b),xmax=b+se_b,alpha=qval<0.01,col=subset))+
         geom_vline(xintercept=0)+
         geom_hline(yintercept=1e-25)+
         geom_text(nudge_x = .3,cex=2.3,data=unique(WTannot[WTannot$set=="Chaperone"&(!is.na(WTannot$set))&(WTannot$subset!=""),-19]),aes(b,-log10(qval),label=symbol))+
        scale_color_brewer(palette = "Dark2")+xlim(-4,NA)
)


ggsave(file="~/GitHub/RNAseq_pipe/Ubiq.pdf",height=3,width=4,
ggplot(WTannot[WTannot$set=="Ubiquitin"&(!is.na(WTannot$set)),])+
  geom_pointrangeh(aes(b,-log10(qval),xmin=(b-se_b),xmax=b+se_b,alpha=qval<0.01,col=set))+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=1e-25)+
  geom_text(data=unique(WTannot[!is.na(WTannot$set)&WTannot$set=="Ubiquitin"&WTannot$qval<0.01,-19]),aes(b,-log10(qval),label=symbol))        -
)



################



mergedDF<-merge.data.frame(WTannot,annot.DF)
genesDF<-data.frame(SO$obs_norm_filt)
annotDF<-merge(genesDF,annotables::grcm38,by.x = "target_id",by.y="ensgene")
annot.DF<-merge(annotDF,allNames,by="symbol",all.x = T)

annot.DF$age=strsplit2(annot.DF$sample,split = "-")[,1]
annot.DF$diff=strsplit2(annot.DF$sample,split = "-")[,3]
annot.DF$sample=strsplit2(annot.DF$sample,split = "-")[,2]

ggplot(annot.DF[!is.na(annot.DF$set),])+geom_boxplot(aes(set,tpm,col=diff))+scale_y_log10()

ggplot(annot.DF[!is.na(annot.DF$subset),])+geom_boxplot(aes(subset,tpm,col=diff))+scale_y_log10()

#ggsave(width= 6,height=6,filename = "SmoothedRPB.pdf",
ggsave("~/GitHub/RNAseq_pipe/RPB_ChangePlot.pdf",width=7, height=5,
  ggplot(mergedDF[!is.na(mergedDF$subset),])+
    geom_line(aes(diff,scaled_reads_per_base,group=target_id),alpha=0.3)+
    geom_line(data = mergedDF[!is.na(mergedDF$subset)&mergedDF$qval<0.01,],aes(diff,scaled_reads_per_base,group=target_id,col=subset))+
    #geom_violin(mapping = aes(x = factor(diff),y=scaled_reads_per_base,color=subset,group=diff))+
    geom_boxplot(cex=.5,mapping = aes(x = factor(diff),y=scaled_reads_per_base,fill=subset,group=diff),width=0.3)+
    facet_grid(set~subset,scales = "free_y")+
    #geom_smooth(mapping = aes(diff,scaled_reads_per_base,color=subset,group=subset),method = "lm")+
    #coord_cartesian(ylim = c(1,100000))+
    scale_y_log10()+
    theme_bw()+scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")
  )

ggsave("~/GitHub/RNAseq_pipe/TPM_ChangePlot.pdf",width=7, height=5,
       ggplot(mergedDF[!is.na(mergedDF$subset),])+
         geom_line(aes(diff,tpm,group=target_id),alpha=0.3)+
         geom_line(data = mergedDF[!is.na(mergedDF$subset)&mergedDF$qval<0.01,],aes(diff,tpm,group=target_id,col=subset))+
         #geom_violin(mapping = aes(x = factor(diff),y=scaled_reads_per_base,color=subset,group=diff))+
         geom_boxplot(cex=.5,mapping = aes(x = factor(diff),y=tpm,fill=subset,group=diff),width=0.3)+
         facet_grid(set~subset,scales = "free_y")+
         #geom_smooth(mapping = aes(diff,scaled_reads_per_base,color=subset,group=subset),method = "lm")+
         #coord_cartesian(ylim = c(1,100000))+
         scale_y_log10()+
         theme_bw()+scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")
)

means<-dcast(annot.DF,symbol~diff,value.var=c("tpm"),fun.aggregate=mean)
sds<-dcast(annot.DF,symbol~diff,value.var=c("tpm"),fun.aggregate=sd)
annot<-dcast(annot.DF,symbol~diff,value.var=c("subset"),fun.aggregate=function(X){as.character(X[1])})
qvals<-dcast(mergedDF$Z,symbol~diff,value.var=c("qval"),fun.aggregate=mean)

ggsave("~/GitHub/RNAseq_pipe/Diff_meanTPM_Scatter.pdf",width=6,height=5,
ggplot()+geom_abline(slope=1)+theme_bw()+
  geom_point(data=means,aes(Var.2,D),col="black",alpha=0.2)+
  geom_abline(slope=1)+
  geom_point(pch=21,data=means[annot$Var.2!=""&!is.na(as.factor(annot$Var.2)),],cex=2.5,aes(Var.2,D,fill=as.factor(annot$Var.2[annot$Var.2!=""&!is.na(as.factor(annot$Var.2))])))+scale_x_log10()+scale_y_log10()+
  guides(fill=guide_legend(title= "Chaperone Class"))+
  ylab("mean Transcripts/Million, Differentiated mNSCs")+
  xlab("mean Transcripts/Million, Undifferentiated mNSCs")
)

means<-dcast(annot.DF,symbol~diff,value.var=c("scaled_reads_per_base"),fun.aggregate=mean)
sds<-dcast(annot.DF,symbol~diff,value.var=c("scaled_reads_per_base"),fun.aggregate=sd)
annot<-dcast(annot.DF,symbol~diff,value.var=c("subset"),fun.aggregate=function(X){as.character(X[1])})
qvals<-dcast(mergedDF,symbol~diff,value.var=c("qval"),fun.aggregate=mean)

ggsave("~/GitHub/RNAseq_pipe/Diff_meanRPB_Scatter.pdf",width=6,height=5,
       ggplot()+geom_abline(slope=1)+theme_bw()+
         geom_point(data=means,aes(Var.2,D),col="black",alpha=0.2)+
         geom_abline(slope=1)+
         geom_point(pch=21,data=means[annot$Var.2!=""&!is.na(as.factor(annot$Var.2)),],cex=2.5,aes(Var.2,D,fill=as.factor(annot$Var.2[annot$Var.2!=""&!is.na(as.factor(annot$Var.2))])))+scale_x_log10()+scale_y_log10()+
         guides(fill=guide_legend(title= "Chaperone Class"))+
         ylab("mean scaled reads/base, Differentiated mNSCs")+
         xlab("mean scaled reads/base, Undifferentiated mNSCs")
)

ggsave(filename = "~/GitHub/RNAseq_pipe/SmoothedCount.pdf",width=6,height=5,
       ggplot(annot.DF[!is.na(annot.DF$subset),])+
         geom_line(aes(diff,scaled_reads_per_base,group=target_id,col=subset),alpha=0.2)+
         facet_wrap(~subset)+
         geom_smooth(mapping = aes(diff,scaled_reads_per_base,color=subset,group=subset),method = "lm", lwd=2)+scale_y_log10()+theme_pubr()+scale_color_brewer(palette="Spectral")
)

