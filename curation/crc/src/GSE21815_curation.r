
rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE21815_full_pdata.csv",as.is=TRUE,row.names=1)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

#alt_sample_name
curated$alt_sample_name <- uncurated$title

##sample_type
curated$sample_type[1:132] <-"tumor"
curated$sample_type[133:141] <-"adjacentnormal"

#age
tmp <-apply(uncurated,1,getVal,string="age: ")
tmp<-sapply(tmp,"[[",1)
tmp <- sub("age: ","",tmp,fixed=TRUE)
curated$age_at_initial_pathologic_diagnosis <- as.numeric(tmp)

#gender
tmp <- apply(uncurated,1,getVal,string="gender: ")
tmp <- sub("gender: ","",tmp,fixed=TRUE)
tmp[tmp=="M"] <- "m"
tmp[tmp=="F"] <- "f"
curated$gender <- tmp

#ethnicity
curated$ethnicity <- "other"

##T
tmp1<- apply(uncurated,1,getVal,string="tumor-nodes-metastasis (tnm) classification - t: ")
tmp1<- sub("tumor-nodes-metastasis (tnm) classification - t: ","",tmp1,fixed=TRUE)
tmp2<- apply(uncurated,1,getVal,string="tumor-nodes-metastasis classification t: ")
tmp2<- sub("tumor-nodes-metastasis classification t: ","",tmp2,fixed=TRUE)
tmp[(tmp1==1) | (tmp2==1)] <-1
tmp[(tmp1==2) | (tmp2==2)] <-2
tmp[(tmp1==3) | (tmp2==3)] <-3
tmp[(tmp1==4) | (tmp2==4)] <-4
tmp[(tmp1=="Tis") | (tmp2=="Tis")] <- "is"
tmp[(tmp1=="ND") | (tmp2=="ND")] <-NA
curated$T<-tmp

##N
tmp1<- apply(uncurated,1,getVal,string="tumor-nodes-metastasis (tnm) classification - n: ")
tmp1<- sub("tumor-nodes-metastasis (tnm) classification - n: ","",tmp1,fixed=TRUE)
tmp2<- apply(uncurated,1,getVal,string="tumor-nodes-metastasis classification n: ")
tmp2<- sub("tumor-nodes-metastasis classification n: ","",tmp2,fixed=TRUE)
tmp[(tmp1==0) | (tmp2==0)] <-0
tmp[(tmp1==1) | (tmp2==1)] <-1
tmp[(tmp1==2) | (tmp2==2)] <-2
tmp[(tmp1=="ND") | (tmp2=="ND")] <-NA
tmp[5]<-"X"
curated$N<-tmp

##M
tmp1<- apply(uncurated,1,getVal,string="tumor-nodes-metastasis (tnm) classification - m: ")
tmp1<- sub("tumor-nodes-metastasis (tnm) classification - m: ","",tmp1,fixed=TRUE)
tmp2<- apply(uncurated,1,getVal,string="tumor-nodes-metastasis classification m: ")
tmp2<- sub("tumor-nodes-metastasis classification m: ","",tmp2,fixed=TRUE)
tmp[(tmp1==0) | (tmp2==0)] <-0
tmp[(tmp1==1) | (tmp2==1)] <-1
tmp[(tmp1=="ND") | (tmp2=="ND")] <-NA
curated$M<-tmp

##summarystage
tmp1 <-curated$T
tmp2 <-curated$N
tmp3 <-curated$M
tmp[(tmp1<4) & (tmp2==0) & (tmp3==0)]<-"early"
tmp[(tmp1==4)]<-"late"
tmp[(tmp1=="is") & (tmp2==0) & (tmp3==0)]<-"early"
tmp[(tmp2>0) | (tmp3>0)] <-"late"
tmp[5]<-NA
tmp[133:141]<-NA
curated$summarystage <-tmp

#stageall
tmp1 <-curated$T
tmp2 <-curated$N
tmp3 <-curated$M
tmp[(tmp1=="is") & (tmp2==0) & (tmp3==0)]<-0
tmp[((tmp1==1) | (tmp1==2)) & (tmp2==0) & (tmp3==0)]<-1
tmp[((tmp1==3) | (tmp1==4)) & (tmp2==0) & (tmp3==0)]<-2
tmp[((tmp2==1) | (tmp2==2)) & (tmp3==0)]<-3
tmp[tmp3==1]<-4
tmp[5]<-NA
curated$stageall <-tmp
 
#Dstage : a better way is to derive from stageall
# tmp <- apply(uncurated,1,getVal,string="dukes: ")
# tmp <- sub("dukes: ","",tmp,fixed=TRUE)
# tmp[tmp=="ND"]<-NA
tmp[tmp==0]<-NA
tmp[tmp==1]<- NA
tmp[tmp==2]<-"B"
tmp[tmp==3]<-"C"
tmp[tmp==4]<-"D"
curated$Dstage <- tmp
curated$Dstage[4]<-"A"
curated$Dstage[14]<-"B"
curated$Dstage[15]<-"B"
curated$Dstage[16]<-"B"
curated$Dstage[43]<-"B"
curated$Dstage[44]<-"B"
curated$Dstage[54]<-"B"
curated$Dstage[61]<-"A"
curated$Dstage[64]<-"B"
curated$Dstage[69]<-"B"
curated$Dstage[75]<-"B"
curated$Dstage[83]<-"B"
curated$Dstage[91]<-"B"
curated$Dstage[95]<-"A"
curated$Dstage[104]<-"B"
curated$Dstage[117]<-"A"
curated$Dstage[127]<-"B"
curated$Dstage[130]<-"B"
curated$Dstage[131]<-"A"
curated$Dstage[132]<-"B"

##preop_drug_treatment
curated$preop_drug_treatment <- "n"

##G
tmp <- apply(uncurated,1,getVal,string="his type: ")
tmp <- sub("his type: ","",tmp,fixed=TRUE)
tmp[tmp=="ND"]<-NA
tmp[tmp=="wel"]<-"1"
tmp[tmp=="mod"]<-"2"
tmp[tmp=="por"]<-"3"
curated$G <- tmp

##summarygrade
tmp[tmp=="1"]<-"low"
tmp[tmp=="2"]<-"low"
tmp[tmp=="3"]<-"high"
curated$summarygrade <-tmp
curated <- postProcess(curated, uncurated)
curated<-updatedfs(curated)

write.table(curated, row.names=FALSE, file="../curated/GSE21815_curated_pdata.txt",sep="\t")

