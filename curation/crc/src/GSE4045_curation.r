rm(list=ls())
source("../../functions.R")

##uncurated <- read.csv("../uncurated/pythonscript1_output_GSE_14333.csv",as.is=TRUE,row.names=1)
uncurated.raw <- read.csv("../uncurated/GSE4045_full_pdata.csv",as.is=TRUE,row.names=1)

##all data is in the column characteristics_ch1, so parse this into its own table:
uncurated <- strsplit(uncurated.raw$description,split=",")
#uncurated <- strsplit(uncurated.raw$description,split=". ")
uncurated <- do.call(rbind,uncurated)
rownames(uncurated) <- rownames(uncurated.raw)
colnames(uncurated) <- 1:ncol(uncurated)

uncurated <- data.frame(uncurated,stringsAsFactors=FALSE)


##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##--------------------
##start the mappings
##--------------------
##alt_sample_name
##gender
##stageall
##MSS
#MSI
##location
##family_history
##summarygrade
##G

##title -> alt_sample_name
tmp <- uncurated.raw$title
curated$alt_sample_name <- tmp

##X6-> gender
tmp <- uncurated$X6
tmp <- sub("male","m",tmp,fixed=TRUE)
tmp <- sub("fem","f",tmp,fixed=TRUE)
tmp <- sub(" m","m",tmp,fixed=TRUE)
tmp <- sub(" f","f",tmp,fixed=TRUE)
curated$gender <- tmp 


##Dukes becomes stageall 
tmp<-uncurated$X3
tmp <- sub(" Dukes Stage c","3",tmp,fixed=TRUE)
tmp <- sub(" Dukes Stage b","2",tmp,fixed=TRUE)
tmp <- sub(" Dukes Stage d","4",tmp,fixed=TRUE)
curated$stageall <- tmp 

##MSS
tmp<-uncurated$X4
tmp[tmp==" MSS"]<-"y"
tmp[tmp==" MSI"]<-"n"
curated$mss<-tmp

#MSI
tmp<-uncurated$X4
tmp[tmp==" MSS"]<-"n"
tmp[tmp==" MSI"]<-"y"
curated$msi<-tmp

##family_history
tmp <- uncurated$X5
tmp <- sub(" no cancer in the family","n",tmp,fixed=TRUE)
tmp <- sub(" cancer in the family","y",tmp,fixed=TRUE)
tmp[tmp==""] <-NA
tmp[tmp==" "] <-NA
curated$family_history <- tmp

## location
tmp <- uncurated$X7
tmp <- sub(" Proximal Location ","proximal",tmp,fixed=TRUE)
tmp <- sub(" Distal Location ","distal",tmp,fixed=TRUE)
curated$location <- tmp

##summarygrade
tmp <- uncurated$X8
tmp <- sub("Tumor Grade 1","low",tmp,fixed=TRUE)
tmp <- sub("Tumor Grade 2","low",tmp,fixed=TRUE)
tmp <- sub("Tumor Grade 3","high",tmp,fixed=TRUE)
tmp <- sub("conventional colorectal tumor",NA,tmp,fixed=TRUE)
tmp <- sub(" low","low",tmp,fixed=TRUE)
tmp <- sub(" high","low",tmp,fixed=TRUE)
curated$summarygrade <- tmp

##G
tmp <- uncurated$X8
tmp <- sub("Tumor Grade 1","1",tmp,fixed=TRUE)
tmp <- sub("Tumor Grade 2","2",tmp,fixed=TRUE)
tmp <- sub("Tumor Grade 3","3",tmp,fixed=TRUE)
tmp <- sub("conventional colorectal tumor",NA,tmp,fixed=TRUE)
curated$G <- tmp

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE4045_curated_pdata.txt",sep="\t")

