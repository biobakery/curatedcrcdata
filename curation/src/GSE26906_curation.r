rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE26906_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##Questions for Levi:
##1) our template gives us [01] for metastasis.  They have 0 or the actual location where it metastasized to-- do we want this location?
##2) for apc_mutation we have [yn]-- they have two of them (first and second) and values associated with them.  do we want these values?
##3) did I do the mutation_apc column correctly?

##alt_sample_name
##T
##summarystage
##age
##gender
##summarylocation
##M
##mutation_apc


#alt_sample_name
tmp <- uncurated$title
curated$alt_sample_name <-tmp

#T  are all "II"
curated$T <- "2"

curated$sample_type <- "tumor"

##summarystage
curated$summarystage <- "early"

#age
tmp <- uncurated$characteristics_ch1.1
tmp <- gsub("[^\\d]","",tmp,perl=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp

##gender
tmp <- uncurated$characteristics_ch1.2
tmp <- sub("gender: ","",tmp,fixed=TRUE)
tmp[tmp=="M"] <- "m"
tmp[tmp=="F"] <- "f"
curated$gender <- tmp

##summarylocation
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("localisation: ","",tmp,fixed=TRUE)
tmp[tmp=="Right"] <- "r"
tmp[tmp=="Left"] <- "l"
curated$summarylocation <- tmp

##M
tmp <- uncurated$characteristics_ch1.4
tmp <- sub("metastasis: ","",tmp,fixed=TRUE)
tmp[tmp=="LIVER"] <- "1"
tmp[tmp=="LUNG"] <- "1"
tmp[tmp=="LIVER/BONE"] <- "1"
tmp[tmp=="CNS"] <- "1"
tmp[tmp=="BONE"] <- "1"
tmp[tmp=="LIVER/LUNG"] <- "1"
curated$M <- tmp

#mutation_apc
tmp <- uncurated$characteristics_ch1.5
tmp[tmp=="first apc mutation: 0"] <- "0"
tmp <- ifelse(tmp =="0","n","y")
curated$mutation_apc <- tmp


#tmp2 <- edit(curated)
write.table(curated, row.names=FALSE, file="../curated/GSE26906_curated_pdata.txt",sep="\t")

