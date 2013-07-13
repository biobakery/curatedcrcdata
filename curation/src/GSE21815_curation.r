rm(list=ls())

source("functions.R")

uncurated <- read.csv("../uncurated/GSE21815_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

#alt_sample_name
curated$alt_sample_name <- uncurated$title

##sample_type
tmp<-apply(uncurated,1,getVal,string="sample type: ")
tmp<-sub("sample type: ","",tmp)
tmp[tmp=="normal"]<-"adjacentnormal"
curated$sample_type<- tmp

#age
tmp <-apply(uncurated,1,getVal,string="age: ")
tmp<-sapply(tmp,"[[",1)
tmp <- sub("age: ","",tmp,fixed=TRUE)
curated$age_at_initial_pathologic_diagnosis <- as.numeric(tmp)

#gender
tmp <- apply(uncurated,1,getVal,string="gender: ")
tmp <- sub("gender: ","",tmp,fixed=TRUE)
tmp[tmp=="M"] <- "m"
tmp[tmp=="F"] <- "f"
curated$gender <- tmp


##T
tmp<- apply(uncurated,1,getVal,string="tumor-nodes-metastasis (tnm) classification - t: ")
tmp<- sub("tumor-nodes-metastasis (tnm) classification - t: ","",tmp,fixed=TRUE)
tmp[tmp=="ND"]<-NA
curated$T<-tmp
##summarystage
tmp<- apply(uncurated,1,getVal,string="tumor-nodes-metastasis (tnm) classification - t: ")
tmp<- sub("tumor-nodes-metastasis (tnm) classification - t: ","",tmp,fixed=TRUE)
tmp[tmp=="ND"]<-NA
tmp[tmp=="1"] <- "early"
tmp[tmp=="2"] <- "early"
tmp[tmp=="3"] <- "late"
tmp[tmp=="4"] <- "late"
curated$summarystage <- tmp

##N
tmp<- apply(uncurated,1,getVal,string="tumor-nodes-metastasis (tnm) classification - n: ")
tmp<- sub("tumor-nodes-metastasis (tnm) classification - n: ","",tmp,fixed=TRUE)
tmp[tmp=="ND"]<-NA
curated$N<-tmp

##M
tmp<- apply(uncurated,1,getVal,string="tumor-nodes-metastasis (tnm) classification - m: ")
tmp<- sub("tumor-nodes-metastasis (tnm) classification - m: ","",tmp,fixed=TRUE)
tmp[tmp=="ND"]<-NA
curated$M<-tmp


#Dukes stage
tmp <- apply(uncurated,1,getVal,string="dukes: ")
tmp <- sub("dukes: ","",tmp,fixed=TRUE)
tmp[tmp=="ND"]<-NA
tmp[tmp=="A"] <- 1
tmp[tmp=="B"] <- 2
tmp[tmp=="C"] <- 3
tmp[tmp=="D"] <- 4
curated$stageall <- tmp


curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE21815_curated_pdata.txt",sep="\t")

###Confirm that normal is adjacentnormal