rm(list=ls())
source("functions.R")

uncurated <- read.csv("../uncurated/GSE12225-GPL3676_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##--------------------
##start the mappings
##--------------------

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title


#T
tmp <- uncurated$characteristics_ch1.3
Tstg <- sub(".*T(\\d)N.", "\\1", tmp)
Tstg[tmp=="Stage: TisN0"] <-NA
curated$T <- Tstg

curated$sample_type <- "tumor"

#N
tmp <- uncurated$characteristics_ch1.3
Nstg <- sub(".*T\\dN(\\d)", "\\1", tmp)
Nstg[tmp=="Stage: TisN0"] <-"0" 
curated$N <- Nstg

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE12225-GPL3676_curated_pdata.txt",sep="\t")
