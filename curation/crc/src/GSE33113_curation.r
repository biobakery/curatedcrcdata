
rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE33113_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE33113/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

##sample_type
tmp<-apply(uncurated,1,getVal,string="tissue: ")
tmp<-sub("tissue: ","",tmp)
tmp[tmp=="primary tumor resection"]<-"tumor"
tmp[tmp=="normal colon mucosa"]<-"adjacentnormal"
curated$sample_type<-tmp

##stageall
tmp<-curated$sample_type
tmp[tmp=="tumor"]<-2
tmp[tmp=="adjacentnormal"]<-NA
curated$stageall<-tmp

##N and M stages
tmp[tmp==2]<-0
curated$N<-tmp
curated$M<-tmp

##Dstage
tmp<-curated$sample_type
tmp[tmp=="tumor"]<-"B"
tmp[tmp=="adjacentnormal"]<-NA
curated$Dstage<-tmp

##age_at_initial_pathologic_diagnosis
tmp<-uncurated$characteristics_ch1.2
tmp<-sub("age at diagnosis: ","",tmp)
tmp<-as.numeric(sub(",",".",tmp))
tmp[23]<-tmp[1]
tmp[91]<-tmp[15]
tmp[92]<-tmp[6]
tmp[95]<-tmp[4]
tmp[96]<-tmp[25]
tmp<-round(tmp)
curated$age_at_initial_pathologic_diagnosis<-tmp

##gender
tmp<-apply(uncurated,1,getVal,string="Sex: ")
tmp<-sub("Sex: ","",tmp)
tmp[23]<-tmp[1]
tmp[91]<-tmp[15]
tmp[92]<-tmp[6]
tmp[95]<-tmp[4]
tmp[96]<-tmp[25]
curated$gender<-tmp

##days_to_recurrence_or_death
tmp<-uncurated$characteristics_ch1.5
tmp<-sub("time to meta or recurrence: ","",tmp)
tmp[23]<-tmp[1]
tmp[91]<-tmp[15]
tmp[92]<-tmp[6]
tmp[95]<-tmp[5]
tmp[96]<-tmp[25]
curated$days_to_recurrence_or_death<-tmp

curated <- postProcess(curated, uncurated)

write.table(curated, row.names=FALSE, file="../curated/GSE33113_curated_pdata.txt",sep="\t")
