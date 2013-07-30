rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE14095_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")


##title -> alt_sample_name
tmp <- uncurated$title
curated$alt_sample_name <- tmp

curated <- postProcess(curated, uncurated)

write.table(curated, row.names=FALSE, file="../curated/GSE14095_curated_pdata.txt",sep="\t")

