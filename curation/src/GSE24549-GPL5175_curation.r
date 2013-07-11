rm(list=ls())
source("functions.R")

uncurated <- read.csv("../uncurated/GSE24549-GPL5175_full_pdata.csv",as.is=TRUE,row.names=1)

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


##stageall 
tmp<-uncurated$characteristics_ch1.1
tmp <- sub("tumor stage: ","",tmp,fixed=TRUE)
tmp[tmp=="1"]<-"early"
tmp[tmp=="NA"] <-NA
tmp[tmp=="2"] <-"early"
tmp[tmp=="3"] <-"late"
tmp[tmp=="4"]<-"late"
curated$summarystage <- tmp 

##MSS
tmp <- ifelse(uncurated$characteristics_ch1.2=="msi-status: MSS","y","n")
curated$mss <- tmp

##MSI
tmp <- uncurated$characteristics_ch1.2
tmp[tmp=="msi-status: MSI-L"] <- "y"
tmp[tmp=="msi-status: MSI-H"] <-"y"
tmp[tmp=="msi-status: MSS"] <-"n"
tmp[tmp=="msi-status: NA"] <-NA
curated$msi <- tmp


curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE24549-GPL5175_curated_pdata.txt",sep="\t")
