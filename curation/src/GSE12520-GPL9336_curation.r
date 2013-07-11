rm(list=ls())

source("functions.R")

uncurated <- read.csv("../uncurated/GSE12520-GPL9336_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

#alt_sample_name
##gender-- ASK LW: characteristics_ch1 says one gender, then, on same row, characteristics_ch2 says the opposite...which one is it?
##age
##location
##summarylocation??
##msi
##mss
##DukesStage-- ASK LW!
##recurrence_status
##days_to_death
##vital_status



#alt_sample_name
curated$alt_sample_name <- uncurated$title

#age
tmp <- uncurated$characteristics_ch1.1
tmp <- sub("age: ","",tmp,fixed=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp

##gender
tmp <- uncurated$characteristics_ch1
tmp <- sub("gender: ","",tmp,fixed=TRUE)
tmp[tmp=="Male"] <- "m"
tmp[tmp=="Female"] <- "f"
curated$gender <- tmp

curated$sample_type <- "tumor"

##location
tmp <- uncurated$characteristics_ch1.2
tmp <- sub("location: ","",tmp,fixed=TRUE)
tmp[tmp=="Ascending Colon"] <- "ascending"
tmp[tmp=="Ascending colon"] <- "ascending"
tmp[tmp=="Caecum"] <- "caecum"
tmp[tmp=="Sigmoid"] <- "sigmoid"
tmp[tmp=="Rectum"] <- "rectum"
tmp[tmp=="Descending Colon"] <- "descending"
tmp[tmp=="n/a"] <- NA
tmp[tmp=="Transverse"] <- "transverse"
tmp[tmp=="Rectosigmoid"] <- "rectosigmoid"
tmp[tmp=="Hepatic flexure"] <- "hepaticflexure"
curated$location <- tmp

#msi
tmp <-uncurated$characteristics_ch1.3
tmp[tmp=="msi: MSI: MSI-L"] <- "y"
tmp[tmp=="msi: MSI: MSS"] <- "n"
curated$msi <- tmp

##mss
tmp <-uncurated$characteristics_ch1.3
tmp[tmp=="msi: MSI: MSI-L"] <- "n"
tmp[tmp=="msi: MSI: MSS"] <- "y"
curated$mss <- tmp

##DukesStage
tmp <- uncurated$characteristics_ch1.4
tmp <- sub("stage: ","",tmp,fixed=TRUE)
tmp[tmp=="Dukes' A"] <- "1"
tmp[tmp=="Dukes' B"] <- "2"
tmp[tmp=="Dukes' C"] <- "3"
tmp[tmp=="Dukes' D"] <- "4"
tmp[tmp=="n/a"] <- NA
curated$stageall <- tmp

##recurrence_status
tmp <- uncurated$characteristics_ch1.5
tmp <- sub("recurrence: ","",tmp,fixed=TRUE)
tmp[tmp=="No"] <- "norecurrence"
tmp[tmp=="Yes"] <- "recurrence"
tmp[tmp=="n/a"] <- NA
curated$recurrence_status <- tmp

##days_to_death
tmp <- uncurated$characteristics_ch1.6
tmp <- sub("survival: ","",tmp,fixed=TRUE)
tmp[tmp=="n/a"] <- NA
tmp <- as.numeric(tmp)
tmp <- tmp * 30  #months to days
curated$days_to_death <- tmp

##vital_status
tmp <- uncurated$characteristics_ch1.7
tmp[tmp=="alive/dead at the end of recorded survival: Alive"] <- "living"
tmp[tmp=="alive/dead at the end of recorded survival: Dead"] <- "deceased"
tmp[tmp=="alive/dead at the end of recorded survival: n/a"] <- NA
tmp[tmp=="alive/dead at the end of recorded survival: Sensored"] <- "living"
curated$vital_status <- tmp

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE12520-GPL9336_curated_pdata.txt",sep="\t")

