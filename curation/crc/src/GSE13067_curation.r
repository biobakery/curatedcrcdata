rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE13067_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE13067/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##title -> alt_sample_name
tmp <- uncurated$title
curated$alt_sample_name <- tmp

#sample_type
curated$sample_type <- "tumor"

##msi
tmp <- uncurated$characteristics_ch1
tmp[tmp=="MSI"] <- "y"
tmp[tmp=="MSS"] <- "n"
curated$msi <- tmp

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE13067_curated_pdata.txt",sep="\t")
