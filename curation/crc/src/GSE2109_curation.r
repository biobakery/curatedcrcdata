rm(list=ls())
source("../../functions.R")
getVal <- function(x,string){
  output <- x[grep(string,x,fixed=TRUE)]
  if(length(output)==0) output <- NA
  return(output)
}

uncurated <- read.csv("../uncurated/GSE2109_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE2109/RAW"
##Filter only colon samples based on primary_site
tmp<-apply(uncurated,1,getVal,string="Primary Site: ")
tmp<-sub("Primary Site: ","",tmp)
tmp.1<-which(tmp=="Rectosigmoid")
#tmp.1<-c(tmp.1,which(tmp=="Colon (Extralymphatic Lymphoma)"))
tmp.1<-c(tmp.1,which(tmp=="Colon"))
tmp.1<-c(tmp.1,which(tmp=="Rectum"))
uncurated<-uncurated[names(tmp.1[order(tmp.1)]),]

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##--------------------
##start the mappings
##--------------------
##alt_sample_name
curated$alt_sample_name <- uncurated$title

##sample_type
tmp<-apply(uncurated,1,getVal,string="Histology: ")
tmp<-sub("Histology: ","",tmp)
tmp[grepl( "^Metastatic", tmp )] <- "metastatic"
tmp[!grepl( "^Metastatic", tmp )] <- "tumor"
curated$sample_type<-tmp

##primarysite
tmp<-apply(uncurated,1,getVal,string="Primary Site: ")
tmp<-sub("Primary Site: ","",tmp)
tmp[tmp=="Colon"]<-"co"
tmp[tmp=="Rectosigmoid"]<-"re"
tmp[tmp=="Rectum"]<-"re"
curated$primarysite<-tmp

##age_at_initial_pathologic_diagnosis
tmp<-apply(uncurated,1,getVal,string="Patient Age: ")
tmp<-sub("Patient Age: ","",tmp)
tmp[tmp=="20-30"]<-25
tmp[tmp=="30-39"]<-35
tmp[tmp=="30-40"]<-35
tmp[tmp=="40-50"]<-45
tmp[tmp=="40-49"]<-45
tmp[tmp=="50-59"]<-55
tmp[tmp=="50-60"]<-55
tmp[tmp=="60-69"]<-65
tmp[tmp=="60-70"]<-65
tmp[tmp=="70-79"]<-75
tmp[tmp=="70-80"]<-75
tmp[tmp=="80-89"]<-85
tmp[tmp=="80-90"]<-85
tmp[tmp=="90-100"]<-95
curated$age_at_initial_pathologic_diagnosis<-tmp

##gender
tmp<-apply(uncurated,1,getVal,string="Gender: ")
tmp<-sub("Gender: ","",tmp)
tmp[tmp=="Male"]<-"m"
tmp[tmp=="Female"]<-"f"
curated$gender<-tmp

#ethnicity 
tmp<-apply(uncurated,1,getVal,string="Ethnic Background: ")
tmp<-sub("Ethnic Background: ","",tmp)
tmp[tmp=="Caucasian"]<-"caucasian"
tmp[tmp=="African-American"]<-"black"
tmp[tmp=="Hispanic"]<-"hispanic"
tmp[!((tmp=="caucasian")|(tmp=="black")|(tmp=="hispanic"))]<-"other"
curated$ethnicity<-tmp

##Pathological T: -> T
tmp <- apply(uncurated,1,getVal,string="Pathological T:")
tmp <- sub("Pathological T: ","",tmp,fixed=TRUE)
tmp <- sub("[abc]","",tmp)
#tmp[tmp=="is"]<-NA
curated$T <- tmp

##summarystage
tmp <- apply(uncurated,1,getVal,string="Pathological Stage:")
tmp <- sub("Pathological Stage: ","",tmp,fixed=TRUE)
tmp[(tmp=="is")|(tmp=="0")|(tmp=="1")|(tmp=="2A")] <-"early"
tmp[(tmp=="2B")|(tmp=="3A")|(tmp=="3B")|(tmp=="3C")|(tmp=="4")] <-"late"
tmp[54]<-"late"
tmp[149]<-"early"
tmp[177]<-"late"
tmp[271]<-"late"
tmp[330]<-"late"
tmp[351]<-"late"
tmp[tmp=="Unknown"]<-NA
curated$summarystage <- tmp

##Pathological stage: -> stageall
tmp <- apply(uncurated,1,getVal,string="Pathological Stage:")
tmp <- sub("Pathological Stage: ","",tmp,fixed=TRUE)
tmp <- sub("[ABC]","",tmp)
tmp[tmp=="X"] <- NA
tmp[tmp=="Unknown"] <- NA
curated$stageall <- tmp

##Dstage
tmp <- apply(uncurated,1,getVal,string="Pathological Dukes Stage:")
tmp <- sub("Pathological Dukes Stage: ","",tmp,fixed=TRUE)
curated$Dstage <- tmp

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
curated$N <- tmp

##M
tmp <- apply(uncurated,1,getVal,string="Pathological M:")
tmp <- sub("Pathological M: ","",tmp,fixed=TRUE)
curated$M <- tmp

curated <- postProcess(curated, uncurated)
curated<-updatedfs(curated)

write.table(curated, row.names=FALSE, file="../curated/GSE2109_curated_pdata.txt",sep="\t")

