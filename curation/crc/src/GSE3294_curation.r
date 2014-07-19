
rm(list=ls())

source("../../functions.R")

uncurated1 <- read.csv("../uncurated/GSE3294_full_pdata.csv",as.is=TRUE,row.names=1)
uncurated2 <- read.csv("../uncurated/GSE_3294 addtl data.csv",as.is=TRUE)

uncurated2 <- strsplit(uncurated2$features,split=" ")
uncurated2 <- do.call(rbind,uncurated2)
rownames(uncurated2) <- rownames(uncurated2)
colnames(uncurated2) <- 1:ncol(uncurated2)
uncurated2 <- data.frame(uncurated2,stringsAsFactors=FALSE)

uncurated1$title2 <- substring( uncurated1$title, first=7 )

uncurated <- merge( uncurated1, uncurated2, by.x="title2", by.y="X1")

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

curated$alt_sample_name<-uncurated$title

##sample_type
curated$sample_type<-"tumor"

##primarysite
tmp<-apply(uncurated,1,getVal,string="Tumor ")
tmp<-sub("Tumor ","",tmp)
#tmp<-sub("\t","",tmp)
tmp[tmp=="Rectum Colon Tissue"]<-"re"
tmp[tmp!="re"]<-"co"
# tmp[tmp=="right colon tissue"]<-"co"
# tmp[tmp=="Right Colon Tissue"]<-"co"
# tmp[tmp=="Left Colon Tissue"]<-"co"
curated$primarysite<-tmp

##G
gradetmp<- uncurated$X5
gradetmp[gradetmp=="G"]<-"1"
gradetmp[gradetmp=="M"]<-"2"
gradetmp[gradetmp=="P"]<-"3"
curated$G<-gradetmp

##summarygrade
gradetmp[(gradetmp=="1") | (gradetmp=="2")]<-"low"
gradetmp[gradetmp=="3"]<-"high"
curated$summarygrade<-gradetmp

##summarylocation
tmp<-apply(uncurated,1,getVal,string="Tumor ")
tmp<-sub("Tumor ","",tmp)
tmp<-sub("\t","",tmp)
tmp[tmp=="Rectum Colon Tissue"]<-"l"
tmp[tmp=="right colon tissue"]<-"r"
tmp[tmp=="Right Colon Tissue"]<-"r"
tmp[tmp=="Left Colon Tissue"]<-"l"
curated$summarylocation<-tmp

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
tmp[tmp.t==3]<-"early"
tmp[tmp.t==4]<-"late"
curated$summarystage<-tmp

##Dstage - need to add to template
tmp<-apply(uncurated,1,getVal,string="Dukes ")
tmp<-sub("Dukes ","",tmp)
curated$Dstage<-tmp

#stageall
tmp<-apply(uncurated,1,getVal,string="AJCC ")
tmp<-sub("AJCC ","",tmp)
tmp<-sub("\t","",tmp)
tmp[tmp=="stage 3"]<-3
curated$stageall<-tmp

curated <- postProcess(curated, uncurated)

write.table(curated, row.names=FALSE, file="../curated/GSE3294_curated_pdata.txt",sep="\t")

