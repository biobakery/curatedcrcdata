
rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE18105_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE18105/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##--------------------
##start the mappings
##--------------------

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

##sample_type
tmp <- uncurated$characteristics_ch1.1
tmp[tmp=="tissue: cancer, LCM"] <-"tumor"
tmp[tmp=="tissue: normal, homogenized"] <-"adjacentnormal"
tmp[tmp=="tissue: cancer, homogenized"] <-"tumor"
curated$sample_type <- tmp

##M
tmp <- uncurated$characteristics_ch1
tmp[tmp=="metastasis: metastatic recurrence"] <-"0"
tmp[tmp=="metastasis: metastasis"] <-"1"
tmp[tmp=="metastasis: none"] <-"0"
curated$M <- tmp

##stageall
tmp <- uncurated$characteristics_ch1
tmp[tmp=="metastasis: metastatic recurrence"] <- NA
tmp[tmp=="metastasis: metastasis"] <-"4"
tmp[tmp=="metastasis: none"] <- NA
curated$stageall <- tmp

##summarystage
tmp1 <- uncurated$characteristics_ch1
tmp1[tmp1=="metastasis: metastatic recurrence"] <-"late"
tmp1[tmp1=="metastasis: metastasis"] <-"late"
tmp1[tmp1=="metastasis: none"] <- NA
curated$summarystage <- tmp1

##recurrence_status
tmp <- uncurated$characteristics_ch1
tmp[tmp=="metastasis: metastatic recurrence"] <-"recurrence"
tmp[tmp=="metastasis: metastasis"] <-"norecurrence"
tmp[tmp=="metastasis: none"] <-"norecurrence"
curated$recurrence_status <- tmp

##preop_drug_treatment 
curated$preop_drug_treatment<-"n"

#ethnicity
curated$ethnicity <- "other"

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE18105_curated_pdata.txt",sep="\t")
