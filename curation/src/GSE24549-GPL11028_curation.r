rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE24549-GPL11028_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##--------------------
##start the mappings
##--------------------


##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

##stage -> T
tmp <- uncurated$characteristics_ch1.1
tmp <- sub("tumor stage: ","",tmp,fixed=TRUE)
curated$T <- tmp 

curated$sample_type <- "tumor"

##stageall 
tmp<-uncurated$characteristics_ch1.1
tmp <- sub("tumor stage: ","",tmp,fixed=TRUE)
tmp[tmp=="NA"] <-NA
tmp[tmp=="2"] <-"early"
tmp[tmp=="3"] <-"late"
curated$summarystage <- tmp 

##MSS
tmp <- ifelse(uncurated$characteristics_ch1.2=="msi-status: MSS","y","n")
curated$mss <- tmp

##MSI
tmp <- uncurated$characteristics_ch1.2
##ask Levi: Is there a shorter/more concise way of writing this code?  perhaps w/ an if, ifelse structure?
tmp[tmp=="msi-status: MSI-L"] <- "y"
tmp[tmp=="msi-status: MSI-H"] <-"y"
tmp[tmp=="msi-status: MSS"] <-"n"
tmp[tmp=="msi-status: NA"] <-NA
curated$msi <- tmp


##dfs_status
##living_norecurrence OR deceased_or_recurrence
##0 = no recurrence and alive
##1 = might be living or might be dead
tmp <- uncurated$characteristics_ch1.3
tmp[tmp=="disease_specific_survival_event: 1"] <- "deceased_or_recurrence"
tmp[tmp=="disease_specific_survival_event: 0"] <- "living_norecurrence"
tmp[tmp=="disease_specific_survival_event: NA"] <- NA
curated$dfs_status <- tmp

##days_to_recurrence_or_death
##decimal
tmp <- uncurated$characteristics_ch1.4
tmp <- sub("disease_specific_survival_years: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- round(tmp * 365)  #years to days
tmp[tmp=="disease_specific_survival_event: NA"] <- NA
curated$days_to_recurrence_or_death <- tmp


#tmp2 <- edit(curated) 
write.table(curated, row.names=FALSE, file="../curated/GSE24549-GPL11028_curated_pdata.txt",sep="\t")
