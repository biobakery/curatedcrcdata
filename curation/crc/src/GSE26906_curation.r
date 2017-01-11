rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE26906_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE26906/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##alt_sample_name
##T
##summarystage
##age
##gender
##summarylocation
##M
##mutation_apc

#alt_sample_name
tmp <- uncurated$title
curated$alt_sample_name <-tmp

##sample_type
curated$sample_type <-"tumor"

##primarysite
curated$primarysite<-"co"

#stageall 
curated$stageall <- "2"

#N
curated$N <-0

#M
curated$M <-0

#Dstage
curated$Dstage <-"B"

#age
tmp <- uncurated$characteristics_ch1.1
tmp <- gsub("[^\\d]","",tmp,perl=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp

##gender
tmp <- uncurated$characteristics_ch1.2
tmp <- sub("gender: ","",tmp,fixed=TRUE)
tmp[tmp=="M"] <- "m"
tmp[tmp=="F"] <- "f"
curated$gender <- tmp

##summarylocation
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("localisation: ","",tmp,fixed=TRUE)
tmp[tmp=="Right"] <- "r"
tmp[tmp=="Left"] <- "l"
curated$summarylocation <- tmp

##recurrence_status
tmp <- uncurated$characteristics_ch1.4
tmp <- sub("metastasis: ","",tmp,fixed=TRUE)
tmp[tmp=="LIVER"] <- "recurrence"
tmp[tmp=="LUNG"] <- "recurrence"
tmp[tmp=="LIVER/BONE"] <- "recurrence"
tmp[tmp=="CNS"] <- "recurrence"
tmp[tmp=="BONE"] <- "recurrence"
tmp[tmp=="LIVER/LUNG"] <- "recurrence"
tmp[tmp==0]<-"norecurrence"
curated$recurrence_status <- tmp

#mutation_apc
tmp <- uncurated$characteristics_ch1.5
tmp[tmp=="first apc mutation: 0"] <- "0"
tmp <- ifelse(tmp =="0","n","y")
curated$mutation_apc <- tmp

#msi
curated$msi <- "MSS"

#drug_treatment
curated$drug_treatment <-"n"

##preop_drug_treatment
curated$preop_drug_treatment <-"n"

##chemotherapy
curated$chemotherapy<-"n"

curated <- postProcess(curated, uncurated)
curated<-updatedfs(curated)
write.table(curated, row.names=FALSE, file="../curated/GSE26906_curated_pdata.txt",sep="\t")
