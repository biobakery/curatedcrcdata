rm(list=ls())

source("../../functions.R")

uncurated <- read.csv("../uncurated/GSE12945_full_pdata.csv",as.is=TRUE,row.names=1)
celfile.dir <- "../../../DATA/GSE12945/RAW"
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")
##--------------------
##start the mappings
##--------------------

##alt_sample_name
##location
##lymphnodesremoved
##lymphnodesinvaded
##vital_status
##days_to_death
##gender
##age
##G
##summarygrade
##T
##summarystage
##N
##M


##title -> alt_sample_name
tmp <- uncurated$title
curated$alt_sample_name <- tmp

##characteristics_ch1 -> location
tmp <- uncurated$characteristics_ch1
tmp <- sub("TumorLocalization: ","",tmp,fixed=TRUE)
tmp[tmp=="colon flexure right"] <- "hepaticflexure"
tmp[tmp=="rectum middle third"] <- "rectum"
tmp[tmp=="rectum upper third"] <- "rectum"
tmp[tmp=="rectum transition upper to middle third"] <- "rectum"
tmp[tmp=="sigma"] <- "sigmoid"
tmp[tmp=="colon ascending"] <- "ascending"
tmp[tmp=="rectum middle to lower third"] <- "rectum"
tmp[tmp=="rectum middle to upper third"] <- "rectum"
tmp[tmp=="rectum lower third"] <- "rectum"
tmp[tmp=="colon transversum"] <- "transverse"
tmp[tmp=="rectum rektosigmoidaler transition"] <- "rectum"
tmp[tmp=="colon descending-sigm."] <- "descending"
tmp[tmp=="colon descending"] <- "descending"
tmp[tmp=="colon"] <- "co"
curated$location <- tmp

##summarylocation
tmp<-curated$location
tmp[tmp=="rectum"]<-"l"
tmp[tmp=="descending"]<-"l"
tmp[tmp=="hepaticflexure"]<-"r"
tmp[tmp=="sigmoid"]<-"l"
tmp[tmp=="transverse"]<-"r"
tmp[tmp=="ascending"]<-"r"
tmp[tmp=="caecum"]<-"r"
tmp[tmp=="co"]<-NA
curated$summarylocation <-tmp

curated$sample_type <- "tumor"

##lymphnodesremoved
tmp <- uncurated$characteristics_ch1.1
tmp <- sub("LymphNodesRemoved: ","",tmp,fixed=TRUE)
curated$lymphnodesremoved <- tmp

##lymphnodesinvaded
tmp <- uncurated$characteristics_ch1.2
tmp <- sub("LymphNodesInvaded: ","",tmp,fixed=TRUE)
#tmp<-as.numeric(tmp)
#tmp[tmp<10]<-NA
curated$lymphnodesinvaded <- tmp

##vital_status
tmp <- uncurated$characteristics_ch1.3
tmp <- sub("Death: ","",tmp,fixed=TRUE)
##death == 1; living == 0
tmp[tmp=="0"] <- "living"
tmp[tmp=="1"] <- "deceased"
curated$vital_status <- tmp

##recurrence_status
tmp <- uncurated$characteristics_ch1.4
tmp <- sub("Rezidiv: ","",tmp,fixed=TRUE)
tmp[tmp=="0"] <- "norecurrence"
tmp[tmp=="1"] <- "recurrence"
curated$recurrence_status <- tmp

#dfs_status
tmp <-curated$vital_status
tmp[tmp=="deceased"]<-"deceased_or_recurrence"
tmp[tmp=="living"]<-"living_norecurrence"
tmp[6]<-NA
tmp[15]<-"deceased_or_recurrence"
tmp[28]<-"deceased_or_recurrence"
tmp[40]<-"deceased_or_recurrence"
curated$dfs_status <-tmp

##days_to_death
tmp <- uncurated$characteristics_ch1.5
tmp <- sub("OverallSurvival_months: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- tmp * 30  #months to days
curated$days_to_death <- tmp

##days_to_recurrence_or_death
tmp <- uncurated$characteristics_ch1.6
tmp <- sub("TumorFreeSurvival_months: ","",tmp,fixed=TRUE)
tmp <- as.numeric(tmp)
tmp <- tmp * 30  #months to days
curated$days_to_recurrence_or_death <- tmp

##gender
tmp <- uncurated$characteristics_ch1.7
tmp <- sub("Gender: ","", tmp, fixed=TRUE)
tmp[tmp=="0"] <- "f"
tmp[tmp=="1"] <- "m"
curated$gender <- tmp

##age_at_initial_pathologic_diagnosis
tmp <-uncurated$characteristics_ch1.8
tmp <- sub("AgeAtDiagnosis: ","",tmp,fixed=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp

##G
tmp <-uncurated$characteristics_ch1.9
tmp <- sub("Grading: ","",tmp,fixed=TRUE)
curated$G <- tmp

##summarygrade
tmp <-uncurated$characteristics_ch1.9
tmp <- sub("Grading: ","",tmp,fixed=TRUE)
tmp[tmp=="1"] <- "low"
tmp[tmp=="2"] <- "low"
tmp[tmp=="3"] <- "high"
tmp[tmp=="4"] <- "high"
curated$summarygrade <- tmp

##T
tmp <-uncurated$characteristics_ch1.10
tmp <- sub("pT: ","",tmp,fixed=TRUE)
curated$T <- tmp

##N
tmp <-uncurated$characteristics_ch1.11
tmp <- sub("pN: ","",tmp,fixed=TRUE)
curated$N <- tmp

##M
tmp <-uncurated$characteristics_ch1.12
tmp <- sub("pM: ","",tmp,fixed=TRUE)
#tmp[tmp=="x"] <- NA
tmp[tmp=="x"]<-"X"
curated$M <- tmp

# #Dstage
tmp1 <-curated$T
tmp2 <-curated$N
tmp3 <-curated$M
# tmp[(tmp1==1) & (tmp2==0) & (tmp3==0)]<-"A"
# tmp[(tmp1==2) & (tmp2==0) & (tmp3==0)]<-"B"
# tmp[(tmp1==3) & (tmp2==0) & (tmp3==0)]<-"B"
# tmp[(tmp1==4) & (tmp2==0) & (tmp3==0)]<-"B"
# tmp[(tmp2==1) & (tmp3==0)]<-"C"
# tmp[(tmp2==2) & (tmp3==0)]<-"C"
# tmp[tmp3==1]<-"D"
# curated$Dstage<-tmp

##summarystage
tmp[(tmp1<4) & (tmp2==0) & (tmp3==0)]<-"early"
tmp[(tmp2>0) | (tmp3>0)] <-"late"
curated$summarystage <-tmp

##stageall
tmp<-apply(uncurated,1,getVal,string="UICC: ")
tmp<-sub("UICC: ","",tmp)
tmp[tmp=="I"]<-1
tmp[tmp=="II"]<-2
tmp[tmp=="III"]<-3
tmp[tmp=="IV"]<-4
curated$stageall<-tmp

##primarysite
tmp<-uncurated$characteristics_ch1.14
tmp<-sub("RectumOrColon: ","",tmp)
tmp[tmp=="Colon"]<-"co"
tmp[tmp=="Rectum"]<-"re"
curated$primarysite<-tmp

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE12945_curated_pdata.txt",sep="\t")

##Questions:
## What is Rezidiv?
