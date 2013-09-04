rm(list=ls())

source("../../functions.R")
source("TCGA/TCGA_load_clinical_data.r")


celfile.dir <- "../../../DATA/TCGA-RNASeqV2/RAW"

##map between cel files and patient barcodes
celfiles.map <- read.delim("../uncurated/unc.edu_COAD.IlluminaHiSeq_RNASeqV2.1.3.0.sdrf.txt",as.is=TRUE)
celfiles.map <- celfiles.map[,match(c("Comment..TCGA.Barcode.", "Extract.Name"),colnames(celfiles.map))]
celfiles.map$alt_sample_name <- sub("-[0-9]{2}[A-Z]-[0-9]{2}[A-Z]-[0-9]{4}-[0-9]{2}","",celfiles.map[,1])

summary(rownames(uncurated) %in% celfiles.map$alt_sample_name)  #8 FALSE, 568 TRUE
summary(celfiles.map$alt_sample_name %in% rownames(uncurated))  #14 FALSE, 585 TRUE

##Keep only samples for which we have a mapping to the celfile - lose 26 samples
keep.ids <- intersect(celfiles.map$alt_sample_name,rownames(uncurated))
uncurated <- uncurated[match(keep.ids,rownames(uncurated)),]
uncurated$Extract.Name <- celfiles.map[match(rownames(uncurated), celfiles.map[,3]),1]


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
curated$alt_sample_name <- uncurated$Extract.Name  ##Put Extract.Name for alt_sample_name
curated$sample_name <- curated$alt_sample_name
curated$unique_patient_ID <- uncurated$unique_patient_id

##--------------------
##start the curation
##--------------------
source("TCGA/TCGA_curation_all_platforms.r")

write.table(curated, row.names=FALSE, file="../curated/TCGA-RNASeqV2_curated_pdata.txt",sep="\t")
