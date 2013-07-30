rm(list=ls())
source("../../functions.R")

uncurated.raw <- read.csv("../uncurated/GSE3964_full_pdata.csv",as.is=TRUE,row.names=1)
uncurated<-uncurated.raw[15:29,]
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##--------------------
##start the mappings
##--------------------

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

##sample_type
tmp<-apply(uncurated,1,getVal,string="Tissue: ")
tmp[grep("tumor colon",tmp)]<-"tumor"
tmp[tmp!="tumor"]<-"adjacentnormal"
curated$sample_type<-tmp

curated<-postProcess(curated,uncurated) 
write.table(curated, row.names=FALSE, file="../curated/GSE3964_curated_pdata.txt",sep="\t")
