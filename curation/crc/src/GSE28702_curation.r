rm(list=ls())

source("../../functions.R")

uncurated.raw <- read.csv("../uncurated/GSE28702_full_pdata.csv",as.is=TRUE,row.names=1)
uncurated<-uncurated.raw[-(which(uncurated.raw$characteristics_ch1.1=="location: Liver")),]
uncurated<-uncurated[-(which(uncurated$characteristics_ch1.1=="location: Lung")),]


##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

##sample_type
tmp<-apply(uncurated,1,getVal,string="lesion: ")
tmp<-sub("lesion: ","",tmp)
tmp[tmp=="Primary"]<-"tumor"
tmp[tmp=="Metastasis"]<-"metastatic"
curated$sample_type<-tmp


##gender
tmp<-apply(uncurated,1,getVal,string="gender: ")
tmp<-sub("gender: ","",tmp)
tmp[tmp=="F"]<-"f"
tmp[tmp=="M"]<-"m"
curated$gender<-tmp

#location
tmp<-apply(uncurated,1,getVal,string="location: ")
tmp<-sub("location: ","",tmp)
tmp[tmp=="Rectum"]<-"rectum"
tmp[tmp=="Descending"]<-"descending"
tmp[tmp=="Transverse"]<-"transverse"
tmp[tmp=="cecum"]<-"caecum"
tmp[tmp=="Ascending"]<-"ascending"
tmp[tmp=="Sigmoid"]<-"sigmoid"
tmp[tmp=="Peritoneum"]<-NA

##drug_response
tmp<-apply(uncurated,1,getVal,string="mfolfox6: ")
tmp<-sub("mfolfox6: ","",tmp)
tmp[tmp=="responder"]<-"y"
tmp[tmp=="non-responder"]<-"n"
curated$drug_name<-"mfolfox6"
curated$drug_response<-tmp


curated <- postProcess(curated, uncurated)

write.table(curated, row.names=FALSE, file="../curated/GSE28702_curated_pdata.txt",sep="\t")
##Location = Peritoneum?
##training and test data differentiation?
