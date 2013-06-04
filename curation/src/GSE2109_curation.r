rm(list=ls())
source("../../functions.R")
getVal <- function(x,string){
  output <- x[grep(string,x,fixed=TRUE)]
  if(length(output)==0) output <- NA
  return(output)
}

uncurated <- read.csv("../uncurated/GSE2109_full_pdata.csv",as.is=TRUE,row.names=1)
uncurated <- uncurated[grep("colon",uncurated$source_name_ch1,ignore.case=TRUE),]


##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##--------------------
##start the mappings
##--------------------
##QUESTIONS FOR LEVI:
##1) Consider GSM231959: it has a category called "T" and then it also has something called Pathological stage...and the values are slightly different.  as it stands, i do not take in pathological stage....what does he recomend?  
##2) should i keep substage all ? (currently commented out to pass error check...)



##alt_sample_name
curated$alt_sample_name <- uncurated$title


##Pathological T: -> T
tmp <- apply(uncurated,1,getVal,string="Pathological T:")
tmp <- sub("Pathological T: ","",tmp,fixed=TRUE)
tmp <- sub("[abc]","",tmp)
tmp[tmp=="X"] <- NA
tmp[tmp=="is"] <-NA
curated$T <- tmp

##Pathological stage: -> stageall
tmp <- apply(uncurated,1,getVal,string="Pathological Stage:")
tmp <- sub("Pathological Stage: ","",tmp,fixed=TRUE)
tmp <- sub("[ABC]","",tmp)
tmp[tmp=="X"] <- NA
tmp[tmp=="0"] <-NA
tmp[tmp=="Unknown"] <- NA
curated$stageall <- tmp

##Pathological stage: -> substageall
#tmp <- apply(uncurated,1,getVal,string="Pathological Stage:")
#tmp <- sub("Pathological Stage: ","",tmp,fixed=TRUE)
#tmp <- gsub("[^\\D]","",tmp,perl=TRUE)
#tmp[tmp=="X"] <- NA
#tmp[tmp==""] <- NA
#tmp[tmp=="A"] <-"a"
#tmp[tmp=="B"] <-"b"
#tmp[tmp=="C"] <-"c"
#tmp[tmp=="Unknown"] <- NA
#curated$substageall <- tmp

##Pathological stage: -> summarystage
tmp <- apply(uncurated,1,getVal,string="Pathological Stage:")
tmp <- sub("Pathological Stage: ","",tmp,fixed=TRUE)
tmp <- sub("[ABC]","",tmp)
tmp[tmp=="X"] <-NA
tmp[tmp=="1"] <-"early"
tmp[tmp=="2"] <-"early"
tmp[tmp=="3"] <-"late"
tmp[tmp=="4"] <-"late"
tmp[tmp=="Unknown"] <- NA
tmp[tmp=="0"] <-NA
curated$summarystage <- tmp

##G
tmp <- apply(uncurated,1,getVal,string="Pathological Grade:")
tmp <- sub("Pathological Grade: ","",tmp,fixed=TRUE)
tmp[tmp=="X"] <- NA
curated$G <- tmp

##summarygrade
tmp <- apply(uncurated,1,getVal,string="Pathological Grade:")
tmp <- sub("Pathological Grade: ","",tmp,fixed=TRUE)
tmp[tmp=="X"] <- NA
tmp[tmp=="1"] <-"low"
tmp[tmp=="2"] <-"low"
tmp[tmp=="3"] <-"high"
tmp[tmp=="4"] <-"high"
curated$summarygrade <- tmp
 
##N 
tmp <- apply(uncurated,1,getVal,string="Pathological N:")
tmp <- sub("Pathological N: ","",tmp,fixed=TRUE)
tmp[tmp=="X"] <- NA
curated$N <- tmp

##M
tmp <- apply(uncurated,1,getVal,string="Pathological M:")
tmp <- sub("Pathological M: ","",tmp,fixed=TRUE)
tmp[tmp=="X"] <- NA
curated$M <- tmp

##age_at_initial_pathologic_diagnosis
tmp <- apply(uncurated,1,getVal,string="Patient Age:")
tmp <- sub("Patient Age: ","",tmp,fixed=TRUE)
tmp[tmp=="60-70"] <- "65"
tmp[tmp=="50-60"] <- "55"
tmp[tmp=="20-30"] <- "25"
tmp[tmp=="70-80"] <- "75"
tmp[tmp=="50-59"] <- "55"
tmp[tmp=="40-50"] <- "45"
tmp[tmp=="80-90"] <- "85"
tmp[tmp=="40-49"] <- "45"
tmp[tmp=="30-40"] <- "35"
tmp[tmp=="60-69"] <- "65"
tmp[tmp=="70-79"] <- "75"
tmp[tmp=="80-89"] <- "85"
tmp[tmp=="90-100"] <- "95"
tmp[tmp=="30-39"] <- "35"
curated$age_at_initial_pathologic_diagnosis <- tmp

##family_history
tmp <- apply(uncurated,1,getVal,string="Family History of Cancer?:")
tmp <- sub("Family History of Cancer?: ","",tmp,fixed=TRUE)
tmp[tmp=="Yes"] <- "y"
tmp[tmp=="No"] <- "n"
curated$family_history <- tmp

write.table(curated, row.names=FALSE, file="../curated/GSE2109_curated_pdata.txt",sep="\t")



