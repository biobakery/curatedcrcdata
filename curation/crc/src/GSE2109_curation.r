rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE2109_full_pdata.csv",as.is=TRUE,row.names=1)

##Filter only colon samples based on primary_site
tmp<-apply(uncurated,1,getVal,string="Primary Site: ")
tmp<-sub("Primary Site: ","",tmp)
tmp.1<-which(tmp=="Rectosigmoid")
tmp.1<-c(tmp.1,which(tmp=="Colon (Extralymphatic Lymphoma)"))
tmp.1<-c(tmp.1,which(tmp=="Colon"))
uncurated<-uncurated[names(tmp.1[order(tmp.1)]),]


##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##--------------------
##start the mappings
##--------------------
##alt_sample_name
curated$alt_sample_name <- uncurated$title


##gender
tmp<-apply(uncurated,1,getVal,string="Gender: ")
tmp<-sub("Gender: ","",tmp)
tmp[tmp=="Male"]<-"m"
tmp[tmp=="Female"]<-"f"
curated$gender<-tmp

#ethnicity
tmp<-apply(uncurated,1,getVal,string="Ethnic Background: ")
tmp<-sub("Ethnic Background: ","",tmp)
tmp[tmp=="Caucasian"]<-"caucasian"
tmp[tmp=="African-American"]<-"black"
tmp[tmp=="Hispanic"]<-"hispanic"
tmp[tmp=="American Indian"]<-"other"
tmp[tmp=="Asian"]<-"other"
tmp[tmp=="Hawaiian"]<-"other"
tmp[114]<-"other"
curated$ethnicity<-tmp

##family_history
tmp <- apply(uncurated,1,getVal,string="Family History of Cancer?:")
tmp <- sub("Family History of Cancer?: ","",tmp,fixed=TRUE)
tmp[tmp=="Yes"] <- "y"
tmp[tmp=="No"] <- "n"
curated$family_history <- tmp

##Pathological T: -> T
tmp <- apply(uncurated,1,getVal,string="Pathological T:")
tmp <- sub("Pathological T: ","",tmp,fixed=TRUE)
tmp <- sub("[abc]","",tmp)
tmp[tmp=="X"] <- NA
tmp[tmp=="is"] <-NA
curated$T <- tmp

##Pathological stage: -> stageall
tmp <- apply(uncurated,1,getVal,string="Pathological Stage:")
tmp <- sub("Pathological Stage: ","",tmp,fixed=TRUE)
tmp <- sub("[ABC]","",tmp)
tmp[tmp=="X"] <- NA
tmp[tmp=="0"] <-NA
tmp[tmp=="Unknown"] <- NA
curated$stageall <- tmp

##Pathological stage: -> summarystage
tmp <- apply(uncurated,1,getVal,string="Pathological Stage:")
tmp <- sub("Pathological Stage: ","",tmp,fixed=TRUE)
tmp <- sub("[ABC]","",tmp)
tmp[tmp=="X"] <-NA
tmp[tmp=="1"] <-"early"
tmp[tmp=="2"] <-"early"
tmp[tmp=="3"] <-"late"
tmp[tmp=="4"] <-"late"
tmp[tmp=="Unknown"] <- NA
tmp[tmp=="0"] <-NA
curated$summarystage <- tmp

##G
tmp <- apply(uncurated,1,getVal,string="Pathological Grade:")
tmp <- sub("Pathological Grade: ","",tmp,fixed=TRUE)
tmp[tmp=="X"] <- NA
curated$G <- tmp

##summarygrade
tmp <- apply(uncurated,1,getVal,string="Pathological Grade:")
tmp <- sub("Pathological Grade: ","",tmp,fixed=TRUE)
tmp[tmp=="X"] <- NA
tmp[tmp=="1"] <-"low"
tmp[tmp=="2"] <-"low"
tmp[tmp=="3"] <-"high"
tmp[tmp=="4"] <-"high"
curated$summarygrade <- tmp
 
##N 
tmp <- apply(uncurated,1,getVal,string="Pathological N:")
tmp <- sub("Pathological N: ","",tmp,fixed=TRUE)
tmp[tmp=="X"] <- NA
curated$N <- tmp

##M
tmp <- apply(uncurated,1,getVal,string="Pathological M:")
tmp <- sub("Pathological M: ","",tmp,fixed=TRUE)
tmp[tmp=="X"] <- NA
curated$M <- tmp

curated <- postProcess(curated, uncurated)

write.table(curated, row.names=FALSE, file="../curated/GSE2109_curated_pdata.txt",sep="\t")
##filtered samples by primary_site