
rm(list=ls())

source("../../functions.R")
uncurated <- read.csv("../uncurated/GSE17537_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE17537/RAW"
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

##sample_type
curated$sample_type <- "tumor"

#age
tmp <- apply(uncurated,1,getVal,string="age: ")
tmp <- sub("age: ","",tmp[1,],fixed=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp

#gender
tmp <- apply(uncurated,1,getVal,string="gender: ")
tmp <- sub("gender: ","",tmp,fixed=TRUE)
tmp[tmp=="male"] <- "m"
tmp[tmp=="female"] <- "f"
curated$gender <- tmp

#ethnicity
tmp <- apply(uncurated,1,getVal,string="ethnicity: ")
tmp <- sub("ethnicity: ","",tmp,fixed=TRUE)
tmp[tmp=="other (not caucasian, black, or hispanic)"] <- "other"
curated$ethnicity <- tmp

#ajcc_stage -> stageall
tmp <- apply(uncurated,1,getVal,string="ajcc_stage: ")
tmp <- sub("ajcc_stage: ","",tmp,fixed=TRUE)
curated$stageall <- tmp

#Dstage
tmp[tmp=="1"] <- NA
tmp[tmp=="2"] <- "B"
tmp[tmp=="3"] <- "C"
tmp[tmp=="4"] <- "D"
curated$Dstage <- tmp

#N
tmp <- apply(uncurated,1,getVal,string="ajcc_stage: ")
tmp <- sub("ajcc_stage: ","",tmp,fixed=TRUE)
tmp[tmp=="1"] <- "0"
tmp[tmp=="2"] <- "0"
tmp[tmp=="3"] <- NA
tmp[tmp=="4"] <- NA
curated$N <- tmp

#M
tmp <- apply(uncurated,1,getVal,string="ajcc_stage: ")
tmp <- sub("ajcc_stage: ","",tmp,fixed=TRUE)
tmp[tmp=="1"] <- "0"
tmp[tmp=="2"] <- "0"
tmp[tmp=="3"] <- "0"
tmp[tmp=="4"] <- "1"
curated$M <- tmp

#summarystage
tmp <- apply(uncurated,1,getVal,string="ajcc_stage: ")
tmp <- sub("ajcc_stage: ","",tmp,fixed=TRUE)
tmp[tmp=="1"] <- "early"
tmp[tmp=="2"] <- NA
tmp[tmp=="3"] <- "late"
tmp[tmp=="4"] <- "late"
curated$summarystage <- tmp

#G
tmp <- apply(uncurated,1,getVal,string="grade: ")
tmp <- gsub("[^\\d]","",tmp,perl=TRUE)
tmp[tmp=="21"] <- "2"
curated$G <- tmp

##summarygrade
tmp <- apply(uncurated,1,getVal,string="grade: ")
tmp <- gsub("[^\\d]","",tmp,perl=TRUE)
tmp[tmp=="1"] <- "low"
tmp[tmp=="2"] <- "low"
tmp[tmp=="3"] <- "high"
tmp[tmp=="4"] <- "high"
tmp[tmp=="21"] <- "low"
tmp[tmp==""] <- NA
curated$summarygrade <- tmp

##vital_status
tmp <- apply(uncurated,1,getVal,string="overall_event (death from any cause): ")
tmp <- sub("overall_event (death from any cause): ","",tmp,fixed=TRUE)
tmp[tmp=="no death"] <- "living"
tmp[tmp=="death"] <- "deceased"
curated$vital_status <- tmp

##recurrence_status
tmp <- apply(uncurated,1,getVal,string="dfs_event (disease free survival; cancer recurrence):")
tmp <- sub("dfs_event (disease free survival; cancer recurrence): ","",tmp,fixed=TRUE)
tmp[tmp=="no recurrence"] <- "norecurrence"
curated$recurrence_status <- tmp

##disease_specific_mortality
tmp <- apply(uncurated,1,getVal,string="dss_event (disease specific survival; death from cancer): ")
tmp <- sub("dss_event (disease specific survival; death from cancer): ","",tmp,fixed=TRUE)
tmp[tmp=="no death"] <- "n"
tmp[tmp=="death"] <- "y"
curated$disease_specific_mortality <- tmp

##days_to_death
tmp <- apply(uncurated,1,getVal,string="overall survival follow-up time:")
tmp <- sub("overall survival follow-up time: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- tmp * 30  #months to days
curated$days_to_death <- tmp

##days_to_tumor_recurrence
tmp <- apply(uncurated,1,getVal,string="dfs_time:")
tmp <- sub("dfs_time: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- tmp * 30  #months to days
curated$days_to_tumor_recurrence <- tmp

##primarysite
curated$primarysite <- "co"
curated <- postProcess(curated, uncurated)
curated<-updatedfs(curated)

write.table(curated, row.names=FALSE, file="../curated/GSE17537_curated_pdata.txt",sep="\t")

##Not sure if these apply to this dataset (paper was unclear). Seems from the GEO description that at least primarysite = colon

# ##preop_drug_treatment
# curated$preop_drug_treatment <- "no"






