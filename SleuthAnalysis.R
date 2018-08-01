#Sleuth Analysis
library(reshape2)
library(data.table)
library(biomaRt)
library(biomartr)
library(sleuth)
library(limma)

kalDirs=list.dirs("/Users/ptdolan/Desktop/mNSC_KallistoAnalysis/",full.names = T)

age=strsplit2(kalDirs,"_")[-c(1,14),4]
exp=strsplit2(kalDirs,"_")[-c(1,14),5]
diff<-strsplit2(exp,"-")[,2]

cond<-paste(age,diff)

input<-data.table(sample=paste(age,exp,sep = "-"),age=age,diff=diff,path=kalDirs[c(-1,-14)])

SO<-sleuth_prep(input,num_cores = 2)

SO <- sleuth_fit(SO, ~age+diff, 'full')
SO <- sleuth_fit(SO, ~age, 'age')
SO <-sleuth_fit(SO, ~diff, 'diff')
SO <-sleuth_fit(SO, ~1, 'null')

SO<- sleuth_lrt(SO,null_model =  'age',alt_model =  'full')
SO<- sleuth_lrt(SO,null_model =  'diff',alt_model =  'full')
SO<- sleuth_lrt(SO,null_model =  'null',alt_model =  'full')

Diff_sleuth_table <- sleuth_results(SO, 'age:full', 'lrt', show_all = FALSE)
Age_sleuth_table <- sleuth_results(SO, 'diff:full', 'lrt', show_all = FALSE)
AgeDiff_sleuth_table <- sleuth_results(SO, 'null:full', 'lrt', show_all = FALSE)

Diff_sleuth_significant <- dplyr::filter(Diff_sleuth_table, qval <= 0.05)
Age_sleuth_significant <- dplyr::filter(Age_sleuth_table, qval <= 0.05)
AgeDiff_sleuth_significant <- dplyr::filter(AgeDiff_sleuth_table, qval <= 0.05)

write.csv(Diff_sleuth_significant,file = "mNSC_Diff_Sig-q05.csv")


########################################################################

