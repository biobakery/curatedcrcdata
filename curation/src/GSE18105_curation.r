rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE18105_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##--------------------
##start the mappings
##--------------------


##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

##sample_type
##is there a better way of doing this?
##ex., search for "cancer", if find, then write "tumor", else, write "adjacentnormal"
tmp <- uncurated$characteristics_ch1.1
tmp[tmp=="tissue: cancer, LCM"] <-"tumor"
tmp[tmp=="tissue: normal, homogenized"] <-"adjacentnormal"
tmp[tmp=="tissue: cancer, homogenized"] <-"tumor"
curated$sample_type <- tmp

##M
tmp <- uncurated$characteristics_ch1
tmp[tmp=="metastasis: metastatic recurrence"] <-"1"
tmp[tmp=="metastasis: metastasis"] <-"1"
tmp[tmp=="metastasis: none"] <-"0"
curated$M <- tmp

#tmp2 <- edit(curated) 
write.table(curated, row.names=FALSE, file="../curated/GSE18105_curated_pdata.txt",sep="\t")
