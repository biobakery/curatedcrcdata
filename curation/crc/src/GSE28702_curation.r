
rm(list=ls())

source("../../functions.R")

uncurated1 <- read.csv("../uncurated/GSE28702_full_pdata.csv",as.is=TRUE,row.names=1)
uncurated2 <- read.csv("../uncurated/GSE28702 addtl data.csv",as.is=TRUE)

tmptitle<- uncurated1$title
tmptitle<-sub("FOLFOX responder training, #","",tmptitle)
tmptitle<-sub("FOLFOX non-responder training, #","",tmptitle)
tmptitle<-sub("FOLFOX responder test, #", "", tmptitle)
tmptitle<-sub("FOLFOX non-responder test, #", "", tmptitle)
uncurated1$title2<-tmptitle

uncurated <- merge( uncurated1, uncurated2, by.x="title2", by.y="Sample.number")

rownames(uncurated)<-uncurated$geo_accession
# celfile.dir <- "../../../DATA/GSE28702/RAW"

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$geo_accession

##sample_type
tmp<-apply(uncurated,1,getVal,string="lesion: ")
tmp<-sub("lesion: ","",tmp)
tmp[tmp=="Primary"]<-"tumor"
tmp[tmp=="Metastasis"]<-"metastatic"
curated$sample_type<-tmp

##G
tmpgrade<-uncurated$Grade
tmpgrade[tmpgrade=="well"]<-"1"
tmpgrade[tmpgrade=="mod"]<-"2"
tmpgrade[tmpgrade=="por"]<-"3"
curated$G<-tmpgrade

##summarygrade
tmpsgrade<-curated$G
tmpsgrade[(tmpsgrade=="1") | (tmpsgrade=="2")]<-"low"
tmpsgrade[tmpsgrade=="3"]<-"high"
curated$summarygrade<-tmpsgrade

##age_at_initial_pathologic_diagnosis
curated$age_at_initial_pathologic_diagnosis <-uncurated$Age

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
tmp[tmp=="Cecum"]<-"caecum"
tmp[tmp=="Ascending"]<-"ascending"
tmp[tmp=="Sigmoid"]<-"sigmoid"
tmp[tmp=="Peritoneum"]<-NA
tmp[tmp=="Liver"]<-NA
tmp[tmp=="Lung"]<-NA
curated$location<-tmp

##summarylocation
tmp<-curated$location
tmp[tmp=="rectum"]<-"l"
tmp[tmp=="descending"]<-"l"
tmp[tmp=="transverse"]<-"r"
tmp[tmp=="caecum"]<-"r"
tmp[tmp=="ascending"]<-"r"
tmp[tmp=="sigmoid"]<-"l"
curated$summarylocation <- tmp

##primarysite
tmp<-curated$location
tmp[tmp=="rectum"]<-"re"
tmp[tmp=="descending"]<-"co"
tmp[tmp=="transverse"]<-"co"
tmp[tmp=="caecum"]<-"co"
tmp[tmp=="ascending"]<-"co"
tmp[tmp=="sigmoid"]<-"co"
curated$primarysite <- tmp

##drug_response
tmp<-apply(uncurated,1,getVal,string="mfolfox6: ")
tmp<-sub("mfolfox6: ","",tmp)
tmp[tmp=="responder"]<-"y"
tmp[tmp=="non-responder"]<-"n"
curated$drug_name<-"mfolfox6"
curated$drug_response<-tmp

##drug_treatment
curated$drug_treatment<-"y"

##preop_drug_treatment
curated$preop_drug_treatment<-"n"

##folfox  
curated$folfox<-"y"

##leucovorin
curated$leucovorin<-"y"

##platin  
curated$platin<-"y"  

##chemotherapy  
curated$chemotherapy<-"y"  


curated <- postProcess(curated, uncurated)
curated<-updatedfs(curated)

write.table(curated, row.names=FALSE, file="../curated/GSE28702_curated_pdata.txt",sep="\t")
##Location = Peritoneum?
##training and test data differentiation?
