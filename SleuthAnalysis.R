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
library(sleuth)
library(reshape2)
library(ggplot2)
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
SO<-sleuth_prep(input,target_mapping = DF,aggregation_column = "ensgene", read_bootstrap_tpm=T,extra_bootstrap_summary=T)

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

#Output tables
Diff_sleuth_table <- sleuth_results(SO, 'age:full', 'lrt', show_all = F)
Age_sleuth_table <- sleuth_results(SO, 'diff:full', 'lrt', show_all = FALSE)

Diff_sleuth_table<-merge(Diff_sleuth_table,annotables::grcm38,by.x = "target_id",by.y="ensgene")
Age_sleuth_table<-merge(Age_sleuth_table,annotables::grcm38,by.x = "target_id",by.y="ensgene")

AgeDiff_sleuth_table <- sleuth_results(SO, 'null:full', 'lrt', show_all = FALSE)
Wald_Diff_sleuth_table <- sleuth_results(SO,test = "diffD",test_type ='wt', show_all = F)
cyto
Diff_sleuth_significant <- dplyr::filter(Diff_sleuth_table, qval <= 0.05)
Age_sleuth_significant <- dplyr::filter(Age_sleuth_table, qval <= 0.05)
AgeDiff_sleuth_significant <- dplyr::filter(AgeDiff_sleuth_table, qval <= 0.05)

write.csv(Diff_sleuth_significant,file = "mNSC_Diff_Sig-q05.csv")
write.csv(AgeDiff_sleuth_significant,file = "mNSC_AgeDiff_Sig-q05.csv")

########################################################################
# PHN analysis 
########################################################################

mergedData<-merge(SO$obs_raw,SO$sample_to_covariates,by="sample")
mergedData$refseq_gene<-strsplit2(mergedData$target_id,split = "\\.")[,1]
mergedData<-unique(merge.data.frame(mergedData,mapping,by="refseq_gene"))

ggplot(mergedData)+geom_boxplot(aes(diff,tpm))

