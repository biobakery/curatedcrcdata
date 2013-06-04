rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE4526_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")


##title -> alt_sample_name
tmp <- uncurated$title
curated$alt_sample_name <- tmp

##stageall
##all Dukes c
curated$stageall <- "3"

##recurrence_status (-) = did not develop; (+) = developed recurrence
tmp <- uncurated$description
tmp[tmp=="Recurrence(+)"] <- "recurrence"
tmp[tmp=="Recurrence(-)"] <- "norecurrence"
curated$recurrence_status <- tmp


#tmp2 <- edit(curated)
write.table(curated, row.names=FALSE, file="../curated/GSE4526_curated_pdata.txt",sep="\t")
