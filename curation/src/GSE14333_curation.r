rm(list=ls())
source("functions.R")

uncurated.raw <- read.csv("../uncurated/GSE14333_full_pdata.csv",as.is=TRUE,row.names=1)

##all data is in the column characteristics_ch1, so parse this into its own table:
uncurated <- strsplit(uncurated.raw$characteristics_ch1,split="; ")
uncurated <- do.call(rbind,uncurated)
rownames(uncurated) <- rownames(uncurated.raw)
colnames(uncurated) <- 1:ncol(uncurated)

uncurated <- data.frame(uncurated,stringsAsFactors=FALSE)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##title -> alt_sample_name
tmp <- uncurated.raw$title
curated$alt_sample_name <- tmp

##SummaryLocation
tmp <- uncurated$X1
tmp<-sub("Location: ","",tmp)
tmp <- sub("Right","r",tmp,fixed=TRUE)
tmp[tmp=="Left"] <- "l"
tmp[tmp=="Rectum"] <- NA
tmp[tmp=="Colon"] <- NA
tmp[tmp==""] <- NA
curated$summarylocation <- tmp

#stageall
tmp <- uncurated$X2
tmp<-sub("DukesStage: ","",tmp)
tmp <- sub("A", "1", tmp, fixed=TRUE)
tmp[tmp=="B"] <- "2"
tmp[tmp=="C"] <- "3"
tmp[tmp=="D"] <- "4"
curated$stageall <- tmp

#summarystage
tmp[tmp=="1"] <- "early"
tmp[tmp=="2"] <- "early"
tmp[tmp=="3"] <- "late"
tmp[tmp=="4"] <- "late"
curated$summarystage <- tmp

curated$sample_type <- "tumor"

#age_at_initial_pathologic_diagnosis
tmp <- uncurated$X3
tmp<-sub("Age_Diag: ","",tmp)
tmp <- round(as.numeric(tmp))
curated$age_at_initial_pathologic_diagnosis <- tmp  

#gender
tmp <- uncurated$X4
tmp<-sub("Gender: ","",tmp)
tmp <- sub("M", "m", tmp, fixed=TRUE)
tmp[tmp=="F"] <- "f"
curated$gender <- tmp

#dfs_status
tmp<- uncurated$X6
tmp<-sub("DFS_Cens: ","",tmp)
tmp <- sub("0", "living_norecurrence", tmp, fixed=TRUE)
tmp[tmp=="1"] <- "deceased_or_recurrence"
tmp[tmp=="NA"] <- NA
curated$dfs_status <- tmp

#days_to_recurrence_or_death
tmp <- uncurated$X5
tmp<-sub("DFS_Time: ","",tmp)
tmp<-as.numeric(tmp)
tmp <- tmp * 30  #months to days
curated$days_to_recurrence_or_death <- tmp

curated <- postProcess(curated, uncurated)
write.table(curated, row.names=FALSE, file="../curated/GSE14333_curated_pdata.txt",sep="\t")
