rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE27544_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

#Questions for Levi:
##1) Do we want anything w/ the controls?  or does alt_sample_name give us enough info?  (ie HLA+ vs HLA-)
##2) Does it matter that under characteristics1.1 in the uncurated, for the ones with MSI, it says high-level microsatellite instability?  I did not pull the words "high-level" to my curation...

##alt_sample_name
curated$alt_sample_name <- uncurated$title

curated$sample_type <- "tumor"

##msi
tmp <- uncurated$description
tmp[tmp=="MSI/HLA+"] <- "y"
tmp[tmp=="MSI/HLA-"] <- "y"
tmp[tmp=="MSS/HLA-"] <- "n"
tmp[tmp=="MSS/HLA+"] <- "n"
curated$msi <- tmp

##mss
tmp <- uncurated$description
tmp[tmp=="MSI/HLA+"] <- "n"
tmp[tmp=="MSI/HLA-"] <- "n"
tmp[tmp=="MSS/HLA-"] <- "y"
tmp[tmp=="MSS/HLA+"] <- "y"
curated$mss <- tmp

#tmp2 <- edit(curated)
write.table(curated, row.names=FALSE, file="../curated/GSE27544_curated_pdata.txt",sep="\t")