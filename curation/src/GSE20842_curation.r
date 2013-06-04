rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE20842_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##Questions for Levi:
##NOTE: wildtype -> "n"; mutated -> "y"
##2) It also has a column saying tumor, do we want to include this for CRC?  

##alt_sample_name
##age
##gender
##KRAS

#alt_sample_name
tmp <- uncurated$characteristics_ch1
tmp <- sub("patient id: ","",tmp,fixed=TRUE)
curated$alt_sample_name <-tmp


#sample_type
tmp <- uncurated$characteristics_ch1.4
tmp[tmp=="tissue: tumor"] <- "tumor"
tmp[tmp=="tissue: mucosa"] <- "adjacentnormal"
curated$sample_type <- tmp

#age
tmp <- uncurated$characteristics_ch1.2
tmp <- sub("age: ","",tmp,fixed=TRUE)
tmp <- round(as.numeric(tmp))
curated$age_at_initial_pathologic_diagnosis <- tmp

##gender
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("gender: ","",tmp,fixed=TRUE)
tmp[tmp=="male"] <- "m"
tmp[tmp=="female"] <- "f"
curated$gender <- tmp

##KRAS
tmp <- uncurated$characteristics_ch1.5
tmp <- sub("genome/variation: ","",tmp,fixed=TRUE)
tmp[tmp=="wild type KRAS"] <- "n"
tmp[tmp=="mutated KRAS"] <- "y"
curated$kras <- tmp

#tmp2 <- edit(curated)
write.table(curated, row.names=FALSE, file="../curated/GSE20842_curated_pdata.txt",sep="\t")

