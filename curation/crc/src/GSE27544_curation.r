
rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE27544_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE27544/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##alt_sample_name
curated$alt_sample_name <- uncurated$title

##sample_type
curated$sample_type <- "tumor"

##msi
tmp <- uncurated$description
tmp[tmp=="MSI/HLA+"] <- "MSI"
tmp[tmp=="MSI/HLA-"] <- "MSI"
tmp[tmp=="MSS/HLA-"] <- "MSS"
tmp[tmp=="MSS/HLA+"] <- "MSS"
curated$msi <- tmp

curated<-postProcess(curated,uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE27544_curated_pdata.txt",sep="\t")
