rm(list=ls())

source("functions.R")

uncurated <- read.csv("../uncurated/GSE45270_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

curated <- postProcess(curated, uncurated)

write.table(curated, row.names=FALSE, file="../curated/GSE45270_curated_pdata.txt",sep="\t")