rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE3294_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

curated$alt_sample_name<-uncurated$title

##gender
tmp<-apply(uncurated,1,getVal,string="Sex ")
tmp<-sub("Sex ","",tmp)
tmp<-sub("\t","",tmp)
tmp[tmp=="M"]<-"m"
tmp[tmp=="F"]<-"f"
curated$gender<-tmp

##age_at_initial_pathologic_diagnosis
tmp<-apply(uncurated,1,getVal,string="Age ")
tmp<-sub("Age ","",tmp)
tmp<-sub("\t","",tmp)
curated$age_at_initial_pathologic_diagnosis<-tmp

##TNM
tmp<-apply(uncurated,1,getVal,string="TNM ")
tmp<-sub("TNM ","",tmp)
tmp<-sub("\t+","",tmp)
tmp.t<-gsub("N[[:digit:]]M[[:digit:]]","",tmp)
tmp.t<-sub("T","",tmp.t)
tmp.m<-gsub('T[[:digit:]]N[[:digit:]]',"",tmp)
tmp.m<-sub("M","",tmp.m)
tmp.n<-gsub('T[[:digit:]]',"",tmp)
tmp.n<-gsub('M[[:digit:]]',"",tmp.n)
tmp.n<-sub("N","",tmp.n)
curated$T<-tmp.t
curated$M<-tmp.m
curated$N<-tmp.n

##summarystage
tmp[tmp.t==1]<-"early"
tmp[tmp.t==2]<-"early"
tmp[tmp.t==3]<-"late"
tmp[tmp.t==4]<-"late"
curated$summarystage<-tmp

#stageall
tmp<-apply(uncurated,1,getVal,string="Dukes ")
tmp<-sub("Dukes ","",tmp)
tmp[tmp=="A"]<-1
tmp[tmp=="B"]<-2
tmp[tmp=="C"]<-3
tmp[tmp=="D"]<-4
curated$stageall<-tmp

curated <- postProcess(curated, uncurated)

write.table(curated, row.names=FALSE, file="../curated/GSE3294_curated_pdata.txt",sep="\t")

##Talk to Dr. Levi about a few reamining entries which come up in the uncurated data file
