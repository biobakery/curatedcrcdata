rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE16125-GPL5175_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE16125-GPL5175/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##--------------------
##start the mappings
##--------------------


##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

#sample_type
curated$sample_type <- "tumor"

##gender
tmp <- uncurated$characteristics_ch1
tmp[tmp=="gender: F"] <-"f"
tmp[tmp=="gender: M"] <-"m"
curated$gender <- tmp

##age_at_initial_pathologic_diagnosis
tmp <- uncurated$characteristics_ch1.1
tmp <- sub("age: ","",tmp,fixed=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp

##stageall
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("stage: ","",tmp,fixed=TRUE)
tmp[tmp=="I"] <-"1"
tmp[tmp=="II"] <-"2"
tmp[tmp=="III"] <-"3"
tmp[tmp=="IV"] <-"4"
tmp[tmp=="NA"] <- NA
curated$stageall <- tmp

##summarystage
tmp[tmp=="1"]<-"early"
tmp[tmp=="2"]<- NA
tmp[tmp=="3"]<-"late"
tmp[tmp=="4"]<-"late"
curated$summarystage<-tmp

#Dstage
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("stage: ","",tmp,fixed=TRUE)
tmp[tmp=="I"] <-NA
tmp[tmp=="II"] <-"B"
tmp[tmp=="III"] <-"C"
tmp[tmp=="IV"] <-"D"
tmp[tmp=="NA"] <- NA
curated$Dstage <- tmp

#M
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("stage: ","",tmp,fixed=TRUE)
tmp[tmp=="I"] <-"0"
tmp[tmp=="II"] <-"0"
tmp[tmp=="III"] <-"0"
tmp[tmp=="IV"] <-"1"
tmp[tmp=="NA"] <- NA
curated$M <- tmp

##vital_status
tmp <- uncurated$characteristics_ch1.4
tmp[tmp=="status: 0"] <-"deceased"
tmp[tmp=="status: 1"] <-"living"
tmp[tmp=="status: NA"] <- NA
curated$vital_status <- tmp

#MSI
curated$msi <- "MSS"

#family_history
curated$family_history <- "n"

#preop_drug_treatment
curated$preop_drug_treatment <- "n"

##kras
tmp <- uncurated$characteristics_ch1.7
tmp[tmp=="kras: 1"] <-"y"
tmp[tmp=="kras: 0"] <-"n"
curated$kras <- tmp

##apc
tmp <- uncurated$characteristics_ch1.6
tmp[tmp=="apc: 1"] <-"y"
tmp[tmp=="apc: 0"] <-"n"
curated$mutation_apc <- tmp

##tp53
# tmp <- uncurated$characteristics_ch1.8
# tmp[tmp=="tp53: 1"] <-"y"
# tmp[tmp=="tp53: 0"] <-"n"
# curated$mutation_apc <- tmp

##days_to_death
tmp <- uncurated$characteristics_ch1.5
tmp <- sub("survival: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- round(tmp * 30)  #months to days
curated$days_to_death <- tmp

##tumor_size
tmp <- uncurated$characteristics_ch1.2
tmp <- sub("dimension: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
curated$tumor_size <- tmp

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE16125-GPL5175_curated_pdata.txt",sep="\t")

