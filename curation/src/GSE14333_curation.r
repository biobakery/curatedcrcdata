rm(list=ls())
source("../../functions.R")
getVal <- function(x,string){
  output <- x[grep(string,x,fixed=TRUE)]
  if(length(output)==0) output <- NA
  return(output)
}

##uncurated <- read.csv("../uncurated/pythonscript1_output_GSE_14333.csv",as.is=TRUE,row.names=1)
uncurated.raw <- read.csv("../uncurated/GSE14333_full_pdata.csv",as.is=TRUE,row.names=1)

##all data is in the column characteristics_ch1, so parse this into its own table:
uncurated <- strsplit(uncurated.raw$characteristics_ch1,split="; ")
uncurated <- do.call(rbind,uncurated)
rownames(uncurated) <- rownames(uncurated.raw)
colnames(uncurated) <- 1:ncol(uncurated)

for (i in 1:ncol(uncurated)){
  before.colvals <- uncurated[,i]
  intermediate.colvals <- strsplit(before.colvals,split=":")
  after.colvals <- do.call(rbind,intermediate.colvals)
  colnames(uncurated)[i] <- unique(after.colvals[,1])
  uncurated[,i] <- sub(" ","",after.colvals[,2],fixed=TRUE)
}

uncurated <- data.frame(uncurated,stringsAsFactors=FALSE)
uncurated$Age_Diag <- as.numeric(uncurated$Age_Diag)
uncurated$DFS_Time <- as.numeric(uncurated$DFS_Time)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##title -> alt_sample_name
tmp <- uncurated.raw$title
curated$alt_sample_name <- tmp

## Location
tmp <- uncurated$Location
tmp <- sub("Right","r",tmp,fixed=TRUE)
#tmp[tmp=="Right"] <- "r"
tmp[tmp=="Left"] <- "l"
tmp[tmp=="Rectum"] <- NA
tmp[tmp=="Colon"] <- NA
tmp[tmp==""] <- NA
curated$summarylocation <- tmp

#stageall
tmp <- uncurated$DukesStage
tmp <- sub("A", "1", tmp, fixed=TRUE)
#tmp[tmp=="A"] <- "i"
tmp[tmp=="B"] <- "2"
tmp[tmp=="C"] <- "3"
tmp[tmp=="D"] <- "4"
curated$stageall <- tmp

#summarystage
tmp <- uncurated$DukesStage
tmp <- sub("A", "early", tmp, fixed=TRUE)
#tmp[tmp=="A"] <- "early"
tmp[tmp=="B"] <- "early"
tmp[tmp=="C"] <- "late"
tmp[tmp=="D"] <- "late"
curated$summarystage <- tmp

curated$sample_type <- "tumor"

curated$subtype <- "other"

#age_at_initial_pathologic_diagnosis
tmp <- uncurated$Age_Diag
tmp <- round(as.numeric(tmp))
curated$age_at_initial_pathologic_diagnosis <- tmp  

#gender
tmp <- uncurated$Gender
tmp <- sub("M", "m", tmp, fixed=TRUE)
#tmp[tmp=="M"] <- "m"
tmp[tmp=="F"] <- "f"
curated$gender <- tmp

#dfs_status
tmp<- uncurated$DFS_Cens
tmp <- sub("0", "living_norecurrence", tmp, fixed=TRUE)
#tmp[tmp=="0"] <- "living_norecurrence"
tmp[tmp=="1"] <- "deceased_or_recurrence"
tmp[tmp=="NA"] <- NA
curated$dfs_status <- tmp

#days_to_recurrence_or_death
tmp <- uncurated$DFS_Time
tmp <- tmp * 30  #months to days
curated$days_to_recurrence_or_death <- tmp

#tmp2 <- edit(curated)
write.table(curated, row.names=FALSE, file="../curated/GSE14333_curated_pdata.txt",sep="\t")