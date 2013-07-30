rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE26682-GPL96_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

##age
tmp<-apply(uncurated,1,getVal,string="age: ")
tmp<-sub("age: ","",tmp)
curated$age_at_initial_pathologic_diagnosis<-tmp

##gender
tmp<-apply(uncurated,1,getVal,string="gender: ")
tmp<-sub("gender: ","",tmp)
tmp[tmp=="Female"]<-"f"
tmp[tmp=="Male"]<-"m"
curated$gender<-tmp

##MSI_status
tmp<-apply(uncurated,1,getVal,string="microsatellite instability (msi) status: ")
tmp[tmp=="microsatellite instability (msi) status: Stable [MSS]"]<-"n"
tmp[tmp=="microsatellite instability (msi) status: Unknown"]<-NA
tmp[tmp=="microsatellite instability (msi) status: High [MSI-H]"]<-"y"
tmp[tmp=="microsatellite instability (msi) status: Low [MSI-L]"]<-"y"
curated$msi<-tmp

##MSS_status
tmp<-apply(uncurated,1,getVal,string="microsatellite instability (msi) status: ")
tmp[tmp=="microsatellite instability (msi) status: Stable [MSS]"]<-"y"
tmp[tmp=="microsatellite instability (msi) status: Unknown"]<-NA
tmp[tmp=="microsatellite instability (msi) status: High [MSI-H]"]<-"n"
tmp[tmp=="microsatellite instability (msi) status: Low [MSI-L]"]<-"n"
curated$mss<-tmp

curated <- postProcess(curated, uncurated)

write.table(curated, row.names=FALSE, file="../curated/GSE26682-GPL96_curated_pdata.txt",sep="\t")


