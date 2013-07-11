rm(list=ls())
source("functions.R")

uncurated.raw <- read.csv("../uncurated/GSE2630_full_pdata.csv",as.is=TRUE,row.names=1)

##all data is in the column description, so parse this into its own table:
uncurated <- strsplit(uncurated.raw$description,split=",")
#uncurated <- strsplit(uncurated.raw$description,split=". ")
uncurated <- do.call(rbind,uncurated)
rownames(uncurated) <- rownames(uncurated.raw)
colnames(uncurated) <- 1:ncol(uncurated)

uncurated <- data.frame(uncurated,stringsAsFactors=FALSE)


##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##--------------------
##start the mappings
##--------------------

##title -> alt_sample_name
tmp <- uncurated.raw$title
curated$alt_sample_name <- tmp

##T 
tmp <-uncurated$X3
tmp <- gsub("[^T2T3]","",tmp,perl=TRUE)
tmp[tmp=="T3"] <- "3"
tmp[tmp=="T2"] <- "2"
curated$T <- tmp


##x2 -> age_at_initial_pathological_diagnosis
tmp <-uncurated$X2 
tmp <- gsub("[^\\d]","",tmp,perl=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp 

##Male.patient-> gender
tmp <- uncurated$X1
tmp[tmp=="Male patient"] <- "m"
tmp[tmp=="Female patient"] <- "f"
curated$gender <- tmp

##summarylocation 
tmp <- uncurated$X5
tmp <- sub(".+(right|left).+","\\1",uncurated$X5)
tmp[tmp=="right"] <- "r"
tmp[tmp=="left"] <- "l"
tmp[(tmp!="l") & (tmp!="r")] <- NA
curated$summarylocation <- tmp

##N
curated$N <-"0"

##metastasis (M)
curated$M <- "0"

##recurrence_status
tmp <- uncurated$X6
tmp <- sub(".+(with recurrence|without recurrence).+","\\1",uncurated$X6)
tmp[tmp=="with recurrence"] <- "recurrence"
tmp[tmp=="without recurrence"] <- "norecurrence"
curated$recurrence_status <- tmp

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE2630_curated_pdata.txt",sep="\t")

