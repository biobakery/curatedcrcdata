rm(list=ls())

source("../../functions.R")

getVal <- function(x,string){
  output <- x[grep(string,x,fixed=TRUE)]
  if(length(output)==0) output <- NA
  return(output)
}


uncurated <- read.csv("../uncurated/GSE17537_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##Questions for Levi

##2) Is "dss_time" becomes "days_to_death"*30 accurate?





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
tmp[tmp=="21"] <- "2"
curated$G <- tmp

##summarygrade
tmp <- uncurated$characteristics_ch1.4
tmp <- gsub("[^\\d]","",tmp,perl=TRUE)
tmp[tmp=="1"] <- "low"
tmp[tmp=="2"] <- "low"
tmp[tmp=="3"] <- "high"
tmp[tmp=="4"] <- "high"
tmp[tmp=="21"] <- "low"
tmp[tmp==""] <- NA
curated$summarygrade <- tmp

curated$sample_type <- "tumor"

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

#tmp2 <- edit(curated)
write.table(curated, row.names=FALSE, file="../curated/GSE17537_curated_pdata.txt",sep="\t")
