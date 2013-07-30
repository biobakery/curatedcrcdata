rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE11237_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##--------------------
##start the mappings
##--------------------


##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

curated$sample_type<-"tumor"

##gender
tmp <- uncurated$characteristics_ch1
tmp[tmp=="Gender: Female"] <-"f"
tmp[tmp=="Gender: Male"] <-"m"
curated$gender <- tmp

##stageall
tmp <- uncurated$characteristics_ch1.1
tmp[tmp=="AJCC Stage: I"] <-"1"
tmp[tmp=="AJCC Stage: II"] <-"2"
tmp[tmp=="AJCC Stage: III"] <-"3"
tmp[tmp=="AJCC Stage: IV"] <-"4"
curated$stageall <- tmp

##T
tmp <- uncurated$characteristics_ch1.2
tmp <- sub("Pathological T: ","",tmp,fixed=TRUE)
curated$T <- tmp

##summarystage
##based on Pathological T value, not on AJCC stageall value...check...okay?!
tmp <- uncurated$characteristics_ch1.2
tmp[tmp=="Pathological T: 1"] <- "early"
tmp[tmp=="Pathological T: 2"] <- "early"
tmp[tmp=="Pathological T: 3"] <- "late"
tmp[tmp=="Pathological T: 4"] <- "late"
curated$summarystage <- tmp

##N
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("Pathological N: ","",tmp,fixed=TRUE)
curated$N <- tmp

##M
tmp <- uncurated$characteristics_ch1.4
tmp <- sub("Pathological M: ","",tmp,fixed=TRUE)
curated$M <- tmp

##grade
#wikipedia: 1 = "well differentiated"; 2 == "moderateley differentiated"; 3 == "poorly ..."; 4 == "undifferentiated"
tmp <- uncurated$characteristics_ch1.5
tmp[tmp=="Grade: moderate"] <- "2"
tmp[tmp=="Grade: poor"] <- "3"
tmp[tmp=="Grade: well"] <- "1"
tmp[tmp=="Grade: moderate to poor"] <- "3"
curated$G <- tmp

##summarygrade
tmp <- uncurated$characteristics_ch1.5
tmp[tmp=="Grade: moderate"] <- "low"
tmp[tmp=="Grade: poor"] <- "high"
tmp[tmp=="Grade: well"] <- "low"
tmp[tmp=="Grade: moderate to poor"] <- "high"
curated$summarygrade <- tmp

##location
tmp <- uncurated$characteristics_ch1.6
tmp[tmp=="Tumor Site: Sigmoid"] <- "sigmoid"
tmp[tmp=="Tumor Site: Cecum"] <- "caecum"
tmp[tmp=="Tumor Site: Rectosigmoid"] <- "rectosigmoid"
tmp[tmp=="Tumor Site: Hepatic Flexure"] <- "hepaticflexure"
tmp[tmp=="Tumor Site: Ascending"] <- "ascending"
tmp[tmp=="Tumor Site: Descending"] <- "descending"
curated$location <- tmp

#drug treatment
tmp<-uncurated$characteristics_ch1.7
tmp[tmp=="Drug Treatment: YES"]<-"y"
tmp[tmp=="Drug Treatment: NO"]<-"n"
curated$drug_treatment<-tmp

#drug_name
tmp[tmp=="y"]<-"celecoxib"
tmp[tmp=="n"]<-NA
curated$drug_name<-tmp

curated <- postProcess(curated, uncurated) 
write.table(curated, row.names=FALSE, file="../curated/GSE11237_curated_pdata.txt",sep="\t")