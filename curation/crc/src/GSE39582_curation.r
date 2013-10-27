rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE39582_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE39582/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

##gender
tmp<-apply(uncurated,1,getVal,string="Sex: ")
tmp<-sub("Sex: ","",tmp)
tmp[tmp=="M"]<-"m"
tmp[tmp=="F"]<-"f"
curated$gender<-tmp

curated$stageall<-2

##age_at_initial_pathologic_diagnosis
tmp<-uncurated$characteristics_ch1.2
tmp<-sub("age.at.diagnosis: ","",tmp)
tmp<-round(as.numeric(tmp))
curated$age_at_initial_pathologic_diagnosis<-tmp

##summarystage
tmp<-apply(uncurated,1,getVal,string="tnm.stage: ")
tmp<-sub("tnm.stage: ","",tmp)
curated$T<-tmp

##location
tmp<-uncurated$characteristics_ch1.4
tmp<-sub("tumor.location: ","",tmp)
curated$location<-tmp

curated <- postProcess(curated, uncurated)

write.table(curated, row.names=FALSE, file="../curated/GSE39582_curated_pdata.txt",sep="\t")


##rfs_even, rfs_delay, normalizedbatch....???
