
rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE17536_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE17536/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##alt_sample_name
##age
##gender
##ethnicity
##ajcc_stage
##grade
##summarygrade
##vital_status
##recurrence_status
##days_to_death

#alt_sample_name
curated$alt_sample_name <- uncurated$title

##primarysite
curated$primarysite <- "co"

##sample_type
curated$sample_type <- "tumor"

##preop_drug_treatment
curated$preop_drug_treatment <- "no"

#age
tmp <- uncurated$characteristics_ch1
tmp <- sub("age: ","",tmp,fixed=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp

#gender
tmp <- uncurated$characteristics_ch1.1
tmp <- sub("gender: ","",tmp,fixed=TRUE)
tmp[tmp=="male"] <- "m"
tmp[tmp=="female"] <- "f"
curated$gender <- tmp

#ethnicity
tmp <- uncurated$characteristics_ch1.2
tmp <- sub("ethnicity: ","",tmp,fixed=TRUE)
tmp[tmp=="other (not caucasian, black, or hispanic)"] <- "other"
curated$ethnicity <- tmp

#ajcc_stage -> stageall
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("ajcc_stage: ","",tmp,fixed=TRUE)
curated$stageall <- tmp

#N
tmp[tmp=="1"] <- "0"
tmp[tmp=="2"] <- "0"
tmp[tmp=="3"] <- NA
tmp[tmp=="4"] <- NA
curated$N <- tmp

#M
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("ajcc_stage: ","",tmp,fixed=TRUE)
tmp[tmp=="1"] <- "0"
tmp[tmp=="2"] <- "0"
tmp[tmp=="3"] <- "0"
tmp[tmp=="4"] <- "1"
curated$M <- tmp

#Dstage
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("ajcc_stage: ","",tmp,fixed=TRUE)
tmp[tmp=="1"] <- NA
tmp[tmp=="2"] <- "B"
tmp[tmp=="3"] <- "C"
tmp[tmp=="4"] <- "D"
curated$Dstage <- tmp

#summarystage
tmp <- sub("ajcc_stage: ","",tmp,fixed=TRUE)
tmp[tmp=="1"] <- "early"
tmp[tmp=="2"] <- NA
tmp[tmp=="3"] <- "late"
tmp[tmp=="4"] <- "late"
curated$summarystage <- tmp

#G
tmp <- uncurated$characteristics_ch1.4
tmp <- gsub("[^\\d]","",tmp,perl=TRUE)
curated$G <- tmp

##summarygrade
tmp <- uncurated$characteristics_ch1.4
tmp <- gsub("[^\\d]","",tmp,perl=TRUE)
tmp[tmp=="1"] <- "low"
tmp[tmp=="2"] <- "low"
tmp[tmp=="3"] <- "high"
tmp[tmp=="4"] <- "high"
curated$summarygrade <- tmp

##vital_status
tmp <- uncurated$characteristics_ch1.5
tmp <- sub("overall_event (death from any cause): ","",tmp,fixed=TRUE)
tmp[tmp=="no death"] <- "living"
tmp[tmp=="death"] <- "deceased"
curated$vital_status <- tmp

##recurrence_status
tmp <- uncurated$characteristics_ch1.7
tmp <- sub("dfs_event (disease free survival; cancer recurrence): ","",tmp,fixed=TRUE)
tmp[tmp=="no recurrence"] <- "norecurrence"
curated$recurrence_status <- tmp

##days_to_death
tmp <- uncurated$characteristics_ch1.8
tmp <- sub("overall survival follow-up time: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- tmp * 30  #months to days
curated$days_to_death <- tmp

##disease_specific_mortality
tmp <- uncurated$characteristics_ch1.6
tmp <- sub("dss_event (disease specific survival; death from cancer): ","",tmp,fixed=TRUE)
tmp[tmp=="no death"] <- "n"
tmp[tmp=="death"] <- "y"
curated$disease_specific_mortality <- tmp

##days_to_tumor_recurrence
tmp <- uncurated$characteristics_ch1.10
tmp <- sub("dfs_time: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- tmp * 30  #months to days
curated$days_to_tumor_recurrence <- tmp

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE17536_curated_pdata.txt",sep="\t")


