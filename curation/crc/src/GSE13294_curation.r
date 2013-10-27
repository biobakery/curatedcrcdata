rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE13294_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE13294/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##title -> alt_sample_name
tmp <- uncurated$title
curated$alt_sample_name <- tmp

##msi
tmp <- uncurated$characteristics_ch1
tmp[tmp=="MSI"] <- "y"
tmp[tmp=="MSS"] <- "n"
curated$msi <- tmp

##mss
tmp <- uncurated$characteristics_ch1
tmp[tmp=="MSS"] <- "y"
tmp[tmp=="MSI"] <- "n"
curated$mss <- tmp

##sample_type
curated$sample_type<-"tumor"

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE13294_curated_pdata.txt",sep="\t")
