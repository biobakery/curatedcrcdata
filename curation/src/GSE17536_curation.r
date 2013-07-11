rm(list=ls())

source("functions.R")

uncurated <- read.csv("../uncurated/GSE17536_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

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


curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE17536_curated_pdata.txt",sep="\t")


##Questions:
##Overall survival == days_to_death?
##is there is diff between overall survival, dfs_time, dss_time?