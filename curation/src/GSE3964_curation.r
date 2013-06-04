rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE3964_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##--------------------
##start the mappings
##--------------------

##Not sure how to handle liver metasistic, 1st and 2nd biopsy.
##Add chemo resistance or sensitive values to experiment template so we can access it if need be

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

##gender
tmp <- uncurated$characteristics_ch1.1
tmp[tmp=="Gender: female;"] <-"f"
tmp[tmp=="Gender: male;"] <-"m"
tmp[tmp==""] <- NA
curated$gender <- tmp

##age_at_initial_pathologic_diagnosis
tmp <- uncurated$characteristics_ch1.2
tmp <-gsub("[^\\d]","",tmp,perl=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp

##stageall
tmp <- uncurated$characteristics_ch1.4
tmp[tmp=="Tumor stage: stade IV (UICC);"] <-"4"
tmp[tmp==""] <- NA
curated$stageall<- tmp

#tmp2 <- edit(curated) 
write.table(curated, row.names=FALSE, file="../curated/GSE3964_curated_pdata.txt",sep="\t")
