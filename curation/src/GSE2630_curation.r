rm(list=ls())
source("../../functions.R")
getVal <- function(x,string){
  output <- x[grep(string,x,fixed=TRUE)]
  if(length(output)==0) output <- NA
  return(output)
}

##uncurated <- read.csv("../uncurated/pythonscript1_output_GSE_14333.csv",as.is=TRUE,row.names=1)
uncurated.raw <- read.csv("../uncurated/GSE2630_full_pdata.csv",as.is=TRUE,row.names=1)

##all data is in the column characteristics_ch1, so parse this into its own table:
uncurated <- strsplit(uncurated.raw$description,split=",")
#uncurated <- strsplit(uncurated.raw$description,split=". ")
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


##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="CRC_Template_May_26_2011.csv")

##--------------------
##start the mappings
##--------------------
##alt_sample_name
##gender
##age_at_initial_pathologic_diagnosis
##stageall
##tumor_size
##summarylocation
##recurrence_status

##QUESTIONS FOR LEVI:

##2) Did I curate "location" correctly here.  
##3) For each patient, it says that they were followed for five years.  IS there anywhere I should record this in the curation?
##4) not sure how helpful "tumor_size" is...? should we keep it?  as it stands, the values for this "tumor_size" dataset are dif. from the ones in our template


##title -> alt_sample_name
tmp <- uncurated.raw$title
curated$alt_sample_name <- tmp

##T 
tmp <-uncurated$X3
tmp <- gsub("[^T2T3]","",tmp,perl=TRUE)
tmp[tmp=="T3"] <- "3"
tmp[tmp=="T2"] <- "2"
curated$T <- tmp


##x2 -> age_at_initial_pathological_diagnosis
tmp <-uncurated$X2 
tmp <- gsub("[^\\d]","",tmp,perl=TRUE)
curated$age_at_initial_pathologic_diagnosis <- tmp 

##Male.patient-> gender
tmp <- uncurated$Male.patient
tmp[tmp=="Male patient"] <- "m"
tmp[tmp=="Female patient"] <- "f"
curated$gender <- tmp

##FIX!!!
##X5 -> location
##tmp <- uncurated$X5
#tmp <- gsub("[^rightleft]","",tmp,perl=TRUE)
#tmp[tmp=="co"] <- "sigmoid"
#curated$location

##summarylocation 
tmp <- uncurated$X5
tmp <- sub(".+(right|left).+","\\1",uncurated$X5)
tmp[tmp=="right"] <- "r"
tmp[tmp=="left"] <- "l"
tmp[(tmp!="l") & (tmp!="r")] <- NA
curated$summarylocation <- tmp

##N
curated$N <-"0"

##metastasis (M)
curated$M <- "0"

##recurrence_status
tmp <- uncurated$X6
tmp <- sub(".+(with recurrence|without recurrence).+","\\1",uncurated$X6)
tmp[tmp=="with recurrence"] <- "recurrence"
tmp[tmp=="without recurrence"] <- "norecurrence"
curated$recurrence_status <- tmp

##tumor_size
#check this: picked 3 and 8, becuase highest integers above middle values (ie 2.5, 7.5)
tmp <- uncurated$X4
tmp <- sub(".+(below|over).+","\\1",uncurated$X4)
tmp <- as.numeric(tmp)
tmp[tmp=="below"] <- "3"
tmp[tmp=="over"] <- "8"
tmp[(tmp!="3") & (tmp!="8")] <- NA
curated$tumor_size <- tmp


#tmp2 <- edit(curated) 
write.table(curated, row.names=FALSE, file="../curated/GSE2630_curated_pdata.txt",sep="\t")
