rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE12945_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")
##--------------------
##start the mappings
##--------------------
##Questions for Levi:
##3) How to interpret "0" and "1" again?  IN regards to gender...
##4) How do I make summarylocation?

##LW: Email them about lymph nodes, otherwise if not hear back, assume 0 is no, more than that is yes.
##ask them about gender.  if don't hear back, don't use it.

##alt_sample_name
##location
##lymphnodesremoved
##lymphnodesinvaded
##vital_status
##days_to_death

##gender
##age
##G
##summarygrade
##T
##summarystage
##N
##M


##title -> alt_sample_name
tmp <- uncurated$title
curated$alt_sample_name <- tmp

##characteristics_ch1 -> location
tmp <- uncurated$characteristics_ch1
tmp <- sub("TumorLocalization: ","",tmp,fixed=TRUE)
tmp[tmp=="colon flexure right"] <- "flexureright"
tmp[tmp=="rectum middle third"] <- "rectum"
tmp[tmp=="rectum upper third"] <- "rectum"
tmp[tmp=="rectum transition upper to middle third"] <- "rectum"
tmp[tmp=="sigma"] <- "sigmoid"
tmp[tmp=="colon ascending"] <- "ascending"
tmp[tmp=="rectum middle to lower third"] <- "rectum"
tmp[tmp=="rectum middle to upper third"] <- "rectum"
tmp[tmp=="rectum lower third"] <- "rectum"
tmp[tmp=="colon transversum"] <- "transverse"
tmp[tmp=="rectum rektosigmoidaler transition"] <- "rectum"
tmp[tmp=="colon descending-sigm."] <- "descending"
tmp[tmp=="colon descending"] <- "descending"
tmp[tmp=="colon"] <- "co"
curated$location <- tmp

curated$sample_type <- "tumor"

##lymphnodesremoved
tmp <- uncurated$characteristics_ch1.1
tmp <- sub("LymphNodesRemoved: ","",tmp,fixed=TRUE)
curated$lymphnodesremoved <- tmp

##lymphnodesremoved
tmp <- uncurated$characteristics_ch1.2
tmp <- sub("LymphNodesInvaded: ","",tmp,fixed=TRUE)
curated$lymphnodesinvaded <- tmp

##vital_status
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("Death: ","",tmp,fixed=TRUE)
##death == 1; living == 0
tmp[tmp=="0"] <- "living"
tmp[tmp=="1"] <- "deceased"
curated$vital_status <- tmp

##days_to_death
tmp <- uncurated$characteristics_ch1.5
tmp <- sub("OverallSurvival_months: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- tmp * 30  #months to days
curated$days_to_death <- tmp

##gender
tmp <- uncurated$characteristics_ch1.7
tmp <- sub("Gender: ","", tmp, fixed=TRUE)
tmp[tmp=="0"] <- "f"
tmp[tmp=="1"] <- "m"
curated$gender <- tmp

##age_at_initial_pathologic_diagnosis
tmp <-uncurated$characteristics_ch1.8
tmp <- sub("AgeAtDiagnosis: ","",tmp,fixed=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp

##G
tmp <-uncurated$characteristics_ch1.9
tmp <- sub("Grading: ","",tmp,fixed=TRUE)
curated$G <- tmp

##summarygrade
tmp <-uncurated$characteristics_ch1.9
tmp <- sub("Grading: ","",tmp,fixed=TRUE)
tmp[tmp=="1"] <- "low"
tmp[tmp=="2"] <- "low"
tmp[tmp=="3"] <- "high"
tmp[tmp=="4"] <- "high"
curated$summarygrade <- tmp

##T
tmp <-uncurated$characteristics_ch1.10
tmp <- sub("pT: ","",tmp,fixed=TRUE)
curated$T <- tmp

##summarystage
tmp <-uncurated$characteristics_ch1.10
tmp <- sub("pT: ","",tmp,fixed=TRUE)
tmp[tmp=="1"] <- "early"
tmp[tmp=="2"] <- "early"
tmp[tmp=="3"] <- "late"
tmp[tmp=="4"] <- "late"
curated$summarystage <- tmp

##recurrence_status
overall_survival <-uncurated$characteristics_ch1.5
overall_survival <- sub("OverallSurvival_months: ","",overall_survival,fixed=TRUE)
overall_survival <- as.numeric(overall_survival)
tumor_free_survival <-uncurated$characteristics_ch1.6
tumor_free_survival <- sub("TumorFreeSurvival_months: ","",tumor_free_survival,fixed=TRUE)
tumor_free_survival <- as.numeric(tumor_free_survival)
vital_status <- uncurated$characteristics_ch1.3
vital_status <- sub("Death: ","",vital_status,fixed=TRUE)
##death == 1; living == 0
vital_status[vital_status=="0"] <- "living"
vital_status[vital_status=="1"] <- "deceased"

if (tumor_free_survival < overall_survival){
tmp3 <- "recurrence"
}else if (vital_status == "living" && tumor_free_survival == overall_survival){
tmp3 <- "norecurrence"
}else if (vital_status == "deceased" && tumor_free_survival == overall_survival){ 
tmp3 <- NA
}else{
tmp3 <- NA
}
curated$recurrence_status <- tmp3

##N
tmp <-uncurated$characteristics_ch1.11
tmp <- sub("pN: ","",tmp,fixed=TRUE)
curated$N <- tmp

##M
tmp <-uncurated$characteristics_ch1.12
tmp <- sub("pM: ","",tmp,fixed=TRUE)
tmp[tmp=="x"] <- NA
curated$M <- tmp

#tmp2 <- edit(curated)
write.table(curated, row.names=FALSE, file="../curated/GSE12945_curated_pdata.txt",sep="\t")