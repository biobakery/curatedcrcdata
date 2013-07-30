##TCGA has both original and follow-up clinical data.  This is
##unusually complicated because an additional form ID was added to the
##follow-up patient IDs, original form data are repeated in the
##follow-up dataset, but some patients are missing int he follow-up
##dataset, some information is missing in the original (like
##days_to_recurrence), and other information is missing in the
##follow-up (like histology).  The approach taken here will be to take
##data from the follow-up where available and from the original
##otherwise, and include a patient even if they are present only in
##the original.
original <- read.delim("../uncurated/clinical_patient_coad.txt",as.is=TRUE,row.names=1,na.strings=c("[Not Available]", "[Not Applicable]", "[Pending]"),head=T)

##stromal data:
clinical.slide <-
read.delim("../uncurated/biospecimen_aliquot_coad.txt", as.is=TRUE,
    na.strings=c("[Not Available]", "[Not Reported]", "[Not Applicable]", "[Pending]"))

# clinical.slide$bcr_sample_barcode <-
#     as.factor(clinical.slide$bcr_sample_barcode)

## TCGA provides stromal for BOTTOM and TOP location, take the mean of these
## to values
# percent_stromal_cells <- sapply(levels(clinical.slide$bcr_sample_barcode),
#     function(id)
#     mean(clinical.slide[clinical.slide$bcr_sample_barcode==id,c("percent_stromal_cells")],na.rm=TRUE))
# 
# percent_normal_cells <- sapply(levels(clinical.slide$bcr_sample_barcode),
#     function(id)
#     mean(clinical.slide[clinical.slide$bcr_sample_barcode==id,c("percent_normal_cells")],na.rm=TRUE))
# 
# percent_tumor_cells <- sapply(levels(clinical.slide$bcr_sample_barcode),
#     function(id)
#     mean(clinical.slide[clinical.slide$bcr_sample_barcode==id,c("percent_tumor_cells")],na.rm=TRUE))

# ## turn NaNs into NAs
# percent_normal_cells[is.na( percent_normal_cells ) ] <- NA    
# percent_tumor_cells[is.na( percent_tumor_cells ) ] <- NA    
# percent_stromal_cells[is.na( percent_stromal_cells ) ] <- NA    


uncurated<-original