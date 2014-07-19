
rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE24550-GPL5175_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE24550-GPL5175/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##--------------------
##start the mappings
##--------------------

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

##stageall
tmp <- uncurated$characteristics_ch1.1
tmp <- sub("tumor stage: ","",tmp,fixed=TRUE)
curated$stageall <- tmp 

##summarystage
tmp<-curated$stageall
tmp[tmp==3]<-"late"
tmp[tmp==2]<-NA
curated$summarystage<-tmp

##Dstage
tmp<-curated$stageall
tmp[tmp==3]<-"C"
tmp[tmp==2]<-"B"
curated$Dstage<-tmp

##drug_treatment
tmp <-curated$stageall
tmp[tmp==2] <-"n"
tmp[tmp==3] <-"y"
curated$drug_treatment <- tmp

#sample_type
tmp <- uncurated$characteristics_ch1
tmp[tmp=="tissue: colorectal cancer biopsy"] <-"tumor"
tmp[tmp=="tissue: normal colonic mucosa biopsy"] <-"adjacentnormal"
curated$sample_type <- tmp

##MSI
tmp <- uncurated$characteristics_ch1.2
tmp[tmp=="msi-status: MSI-L"] <- "MSS"
tmp[tmp=="msi-status: MSI-H"] <-"MSI"
tmp[tmp=="msi-status: MSS"] <-"MSS"
tmp[tmp=="msi-status: NA"] <-NA
curated$msi <- tmp

##preop_drug_treatment
curated$preop_drug_treatment <- "n" 

##M
curated$M <-0
curated$M[78:90]<-NA

##dfs_status - note: in the uncurated dataset, ch1.3 and 1.4 variable names are switched
tmp <- uncurated$characteristics_ch1.4
tmp <- sub("disease_free_survival_years: ","",tmp,fixed=TRUE)
tmp[tmp==0] <-"living_norecurrence"
tmp[tmp==1] <-"deceased_or_recurrence"
curated$dfs_status <- tmp 

##days_to_recurrence_or_death
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("disease_free_survival_event: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- tmp * 365
curated$days_to_recurrence_or_death <-tmp

#chemotherapy
tmp <-curated$stageall
tmp[tmp==2] <-"n"
tmp[tmp==3] <-"y"
curated$chemotherapy <- tmp

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE24550-GPL5175_curated_pdata.txt",sep="\t")





