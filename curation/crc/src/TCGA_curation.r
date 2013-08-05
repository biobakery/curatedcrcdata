rm(list=ls())

source("../../functions.R")
source("TCGA/TCGA_load_clinical_data.r")
celfile.dir <- "../../../DATA/TCGA/RAW"

##map between cel files and patient barcodes
celfiles.map <- read.delim("../uncurated/unc.edu_COAD.AgilentG4502A_07_3.sdrf.txt",as.is=TRUE)
celfiles.map <- celfiles.map[,match(c("Extract.Name","Array.Data.File"),colnames(celfiles.map))]
celfiles.map$alt_sample_name <- celfiles.map[,1]
celfiles.map<-celfiles.map[which(celfiles.map[,"Extract.Name"]!="Stratagene Univeral Reference"),]

##batch information
batch.info <- readLines("../uncurated/TCGA_file_sources.txt")
batch.info <- do.call(rbind,strsplit(batch.info,split="/"))
batch.info <- batch.info[,-1]
batch.info[,1] <- sub("^.+Level_3\\.","",batch.info[,1])
batch.info[,1] <- sub("\\.5\\.0","",batch.info[,1])
batch.info[,2]<-sub("_lmean.out.logratio.gene.tcga_level3.data.txt","",batch.info[,2])
batch.info <- data.frame(batch.info,stringsAsFactors=FALSE)
colnames(batch.info) <- c("batch","Array.Data.File")
batch.info <- batch.info[match(celfiles.map$Array.Data.File,batch.info$Array.Data.File),]

if(identical(all.equal(celfiles.map$Array.Data.File,batch.info$Array.Data.File),TRUE))
{
  print("Adding batch info to celfiles.map")
  celfiles.map$batch <- batch.info$batch
}

summary(clinical.slide$bcr_aliquot_uuid %in% celfiles.map$alt_sample_name)  
summary(celfiles.map$alt_sample_name %in% clinical.slide$bcr_aliquot_uuid)  

##Keep only samples for which we have a mapping to the celfile - lose 26 samples
keep.ids <- intersect(celfiles.map$alt_sample_name,clinical.slide$bcr_aliquot_uuid)
keep.ids<-clinical.slide[which(clinical.slide[,"bcr_aliquot_uuid"] %in% keep.ids),"bcr_sample_barcode"]
Extract.Name<-keep.ids
keep.ids<-substr(keep.ids,1,12)
uncurated <- uncurated[match(keep.ids,rownames(uncurated)),]
uncurated$Extract.Name<-Extract.Name

##Keep primary tumors and normal tissues only.
tumor.num <- strsplit(uncurated$Extract.Name,split="-")
tumor.num <- sapply(tumor.num,function(x) x[4])
tumor.num <- sub("[A-Z]","",tumor.num)
##This excludes some type 02 - "Recurrent Solid Tumor" and 11 - "Solid
##Tissue Normal".  Type 01 is "Primary Solid Tumor".
keep.samples <- tumor.num=="01" | tumor.num=="11"  ##keep 01 - "Recurrent Solid Tumor" and 11 - "Solid Tissue Normal".
uncurated <- uncurated[keep.samples,]
tumor.num <- tumor.num[keep.samples]

uncurated$unique_patient_id <- rownames(uncurated)
rownames(uncurated) <- make.names(uncurated$Extract.Name)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")
curated$batch <- uncurated$batch
curated$alt_sample_name <- uncurated$Extract.Name  ##Put Extract.Name for alt_sample_name
curated$batch <- uncurated$batch
curated$unique_patient_ID <- uncurated$unique_patient_id

##--------------------
##start the curation
##--------------------
source("TCGA/TCGA_curation_all_platforms.r")

write.table(curated, row.names=FALSE, file="../curated/TCGA_curated_pdata.txt",sep="\t")
