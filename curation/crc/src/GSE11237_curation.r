
rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE11237_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE11237/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##--------------------
##start the mappings
##--------------------


##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

curated$sample_type<-"tumor"

##gender
tmp <- uncurated$characteristics_ch1
tmp[tmp=="Gender: Female"] <-"f"
tmp[tmp=="Gender: Male"] <-"m"
curated$gender <- tmp

##stageall
tmp <- uncurated$characteristics_ch1.1
tmp[tmp=="AJCC Stage: I"] <-"1"
tmp[tmp=="AJCC Stage: II"] <-"2"
tmp[tmp=="AJCC Stage: III"] <-"3"
tmp[tmp=="AJCC Stage: IV"] <-"4"
curated$stageall <- tmp

##T
tmp <- uncurated$characteristics_ch1.2
tmp <- sub("Pathological T: ","",tmp,fixed=TRUE)
curated$T <- tmp

##N
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("Pathological N: ","",tmp,fixed=TRUE)
curated$N <- tmp

##M
tmp <- uncurated$characteristics_ch1.4
tmp <- sub("Pathological M: ","",tmp,fixed=TRUE)
curated$M <- tmp

##summarystage
tmp1 <-curated$T
tmp2 <-curated$N
tmp3 <-curated$M
tmp[tmp1==4]<-"late"
tmp[(tmp1<4) & (tmp2==0) & (tmp3==0)]<-"early"
tmp[(tmp2>0) | (tmp3>0)] <-"late"
curated$summarystage <-tmp

##grade
tmp <- uncurated$characteristics_ch1.5
tmp[tmp=="Grade: moderate"] <- "2"
tmp[tmp=="Grade: poor"] <- "3"
tmp[tmp=="Grade: well"] <- "1"
tmp[tmp=="Grade: moderate to poor"] <- "3"
curated$G <- tmp

##summarygrade
tmp <- uncurated$characteristics_ch1.5
tmp[tmp=="Grade: moderate"] <- "low"
tmp[tmp=="Grade: poor"] <- "high"
tmp[tmp=="Grade: well"] <- "low"
tmp[tmp=="Grade: moderate to poor"] <- "high"
curated$summarygrade <- tmp

##location
tmp <- uncurated$characteristics_ch1.6
tmp[tmp=="Tumor Site: Sigmoid"] <- "sigmoid"
tmp[tmp=="Tumor Site: Cecum"] <- "caecum"
tmp[tmp=="Tumor Site: Rectosigmoid"] <- "rectosigmoid"
tmp[tmp=="Tumor Site: Hepatic Flexure"] <- "hepaticflexure"
tmp[tmp=="Tumor Site: Ascending"] <- "ascending"
tmp[tmp=="Tumor Site: Descending"] <- "descending"
curated$location <- tmp

##summarylocation
tmp <- uncurated$characteristics_ch1.6
tmp[tmp=="Tumor Site: Cecum"] <- "r"
tmp[tmp=="Tumor Site: Hepatic Flexure"] <- "r"
tmp[tmp=="Tumor Site: Ascending"] <- "r"
tmp[tmp!="r"]<-"l"
curated$summarylocation <- tmp

##preop drug treatment
tmp<-uncurated$characteristics_ch1.7
tmp[tmp=="Drug Treatment: YES"]<-"y"
tmp[tmp=="Drug Treatment: NO"]<-"n"
curated$preop_drug_treatment<-tmp
 
##preop drug_name
tmp[tmp=="y"]<-"celecoxib"
tmp[tmp=="n"]<-NA
curated$preop_drug_name<-tmp

curated <- postProcess(curated, uncurated) 
write.table(curated, row.names=FALSE, file="../curated/GSE11237_curated_pdata.txt",sep="\t")
