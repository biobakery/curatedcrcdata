rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE4526_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##title -> alt_sample_name
tmp <- uncurated$title
curated$alt_sample_name <- tmp

##sample_type
curated$sample_type <- "tumor"

##stageall
curated$stageall <- "3"

##summarystage
curated$summarystage <- "late"

##M
curated$M<-"0"

##recurrence_status (-) = did not develop; (+) = developed recurrence
tmp <- uncurated$description
tmp[tmp=="Recurrence(+)"] <- "recurrence"
tmp[tmp=="Recurrence(-)"] <- "norecurrence"
curated$recurrence_status <- tmp

##ethnicity
curated$ethnicity <- "other"

##Dstage
curated$Dstage <- "C"

##preop_drug_treatment
curated$preop_drug_treatment <-"n"

curated <- postProcess(curated, uncurated)

write.table(curated, row.names=FALSE, file="../curated/GSE4526_curated_pdata.txt",sep="\t")
