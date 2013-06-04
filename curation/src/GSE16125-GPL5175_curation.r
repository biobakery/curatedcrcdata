rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE16125-GPL5175_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##--------------------
##start the mappings
##--------------------


##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

curated$sample_type <- "tumor"

curated$msi <- "y"

##gender
tmp <- uncurated$characteristics_ch1
tmp[tmp=="gender: F"] <-"f"
tmp[tmp=="gender: M"] <-"m"
curated$gender <- tmp

##age_at_initial_pathologic_diagnosis
tmp <- uncurated$characteristics_ch1.1
tmp <- sub("age: ","",tmp,fixed=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp

##T
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("stage: ","",tmp,fixed=TRUE)
tmp[tmp=="I"] <-"1"
tmp[tmp=="II"] <-"2"
tmp[tmp=="III"] <-"3"
tmp[tmp=="IV"] <-"4"
tmp[tmp=="NA"] <- NA
curated$T <- tmp

##vital_status
tmp <- uncurated$characteristics_ch1.4
#done correctly?!?!
tmp[tmp=="status: 0"] <-"living"
tmp[tmp=="status: 1"] <-"deceased"
tmp[tmp=="status: NA"] <- NA
curated$vital_status <- tmp

##kras
tmp <- uncurated$characteristics_ch1.7
#done correctly?!?!
tmp[tmp=="kras: 1"] <-"y"
tmp[tmp=="kras: 0"] <-"n"
curated$kras <- tmp

##kras
tmp <- uncurated$characteristics_ch1.6
#done correctly?!?!
tmp[tmp=="apc: 1"] <-"y"
tmp[tmp=="apc: 0"] <-"n"
curated$mutation_apc <- tmp

##tp53
tmp <- uncurated$characteristics_ch1.8
#done correctly?!?!
tmp[tmp=="tp53: 1"] <-"y"
tmp[tmp=="tp53: 0"] <-"n"
curated$mutation_apc <- tmp

##days_to_death ?
tmp <- uncurated$characteristics_ch1.5
#done correctly?!?!
tmp <- sub("survival: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- round(tmp * 30)  #months to days
curated$days_to_death <- tmp



#tmp2 <- edit(curated) 
write.table(curated, row.names=FALSE, file="../curated/GSE16125-GPL5175_curated_pdata.txt",sep="\t")
