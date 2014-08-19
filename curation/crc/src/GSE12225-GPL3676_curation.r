
rm(list=ls())
source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE12225-GPL3676_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE12225-GPL3676/RAW"
tmp<-apply(uncurated,1,getVal,string="Group: ")
tmp<-sub("Group: ","",tmp)
tmp.1<-which(!((tmp=="AA")|(tmp=="AC")))
uncurated<-uncurated[names(tmp.1[order(tmp.1)]),]

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##--------------------
##start the mappings
##--------------------

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$title

##primarysite
curated$primarysite <- "re"

#T
tmp <- uncurated$characteristics_ch1.3
Tstg <- sub(".*T(\\d)N.", "\\1", tmp)
Tstg[tmp=="Stage: TisN0"] <- "is"
curated$T <- Tstg

##sample_type
curated$sample_type <- "tumor"

#N
tmp <- uncurated$characteristics_ch1.3
Nstg <- sub(".*T\\dN(\\d)", "\\1", tmp)
Nstg[tmp=="Stage: TisN0"] <-"0"
Nstg[Nstg=="3"]<-"2"
curated$N <- Nstg

##summarystage
tmp <-curated$N
tmp[tmp=="0"]<- "early"
tmp[tmp!="early"]<- "late"
curated$summarystage <- tmp

#M
curated$M <- "0"

#location
curated$location <- "rectum"

#summarylocation
curated$summarylocation <- "l"

#stageall
tmp <- uncurated$characteristics_ch1.3
tmp[tmp=="Stage: T0N0"]<- NA
tmp[tmp=="Stage: TisN0"]<- "0"
tmp[(tmp=="Stage: T1N0")|(tmp=="Stage: T2N0")] <- "1"
tmp[!((tmp=="0")|(tmp=="1"))]<-"3"
curated$stageall<-tmp

# #Dstage
# tmp1 <-curated$T
# tmp2 <-curated$N
# tmp3 <-curated$M
# tmp[(tmp1==1) & (tmp2==0) & (tmp3==0)]<-"A"
# tmp[(tmp1==2) & (tmp2==0) & (tmp3==0)]<-"B"
# tmp[(tmp1==3) & (tmp2==0) & (tmp3==0)]<-"B"
# tmp[(tmp1==4) & (tmp2==0) & (tmp3==0)]<-"B"
# tmp[(tmp2==1) & (tmp3==0)]<-"C"
# tmp[(tmp2==2) & (tmp3==0)]<-"C"
# tmp[tmp3==1]<-"D"
# curated$Dstage<-tmp

#preop_drug_treatment
curated$preop_drug_treatment <- "n"

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE12225-GPL3676_curated_pdata.txt",sep="\t")
