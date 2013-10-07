##01 - "Recurrent Solid Tumor" and 11 - "Solid Tissue Normal".
##histological_type -> sample_type
tmp <- tumor.num
tmp <- sub("01","tumor",tmp,fixed=TRUE)
tmp <- sub("11","adjacentnormal",tmp,fixed=TRUE)
curated$sample_type <- tmp

##primary_site -> tumor_tissue_site
tmp <- uncurated$tumor_tissue_site 
tmp[tmp=="Colon"] <- "co"
tmp[tmp=="Rectum"]<-"re"
tmp[is.na(tmp)] <- NA 
curated$primarysite <- tmp

##age_at_initial_pathologic_diagnosis
tmp <- uncurated$age_at_initial_pathologic_diagnosis
tmp[tmp=="null"] <- NA
tmp <- as.integer(tmp)
curated$age_at_initial_pathologic_diagnosis <- tmp

##days_to_death
daystodeath <- uncurated$days_to_death  
daystolastfollowup <- uncurated$days_to_last_followup
vitalstatus <- uncurated$vital_status
vitalstatus[is.na(vitalstatus)] <- NA   
tmp <- daystodeath
tmp[grep("LIVING", vitalstatus)] <- daystolastfollowup[grep("LIVING", vitalstatus)]
##If vital status is unknown, set days_to_death to NA as well.
tmp[ is.na(vitalstatus) ] <- NA
curated$days_to_death<- tmp

##vital status
tmp<-uncurated$vital_status
tmp[tmp=="LIVING"]<-"living"
tmp[tmp=="DECEASED"]<-"deceased"
tmp[tmp=="Alive"]<-"living"
tmp[tmp=="Dead"]<-"deceased"
curated$vital_status<-tmp

##MSI
tmp<-uncurated$microsatellite_instability
tmp[tmp=="YES"]<-"y"
tmp[tmp=="NO"]<-"n"
curated$msi<-tmp
#MSS
tmp<-uncurated$microsatellite_instability
tmp[tmp=="YES"]<-"n"
tmp[tmp=="NO"]<-"y"
curated$mss<-tmp

##location
tmp<-uncurated$anatomic_neoplasm_subdivision
tmp[tmp=="Sigmoid Colon"]<-"sigmoid"
tmp[tmp=="Transverse Colon"]<-"transverse"
tmp[tmp=="Cecum"]<-"caecum"
tmp[tmp=="Ascending Colon"]<-"ascending"
tmp[tmp=="Hepatic Flexure"]<-"hepaticflexure"
tmp[tmp=="Splenic Flexure"]<-"splenicflexure"
tmp[tmp=="Descending Colon"]<-"descending"
tmp[tmp=="Rectum"]<-"rectum"
curated$location<-tmp

##gender
tmp<-uncurated$gender
tmp[tmp=="MALE"]<-"m"
tmp[tmp=="FEMALE"]<-"f"
curated$gender<-tmp

##kras
tmp<-uncurated$kras_mutation_found
tmp[tmp=="YES"]<-"y"
tmp[tmp=="NO"]<-"n"
curated$kras<-tmp

curated$braf<-NA

##stage
#tmp<-uncurated$pathologic_stage
#tmp[tmp=="Stage IIA"]<-2
#tmp[tmp=="Stage II"]<-2
#tmp[tmp=="Stage IIB"]<-2
#tmp[tmp="Stage IIIB"]<-3
#tmp[tmp="Stage III"]<-3
#tmp[tmp="Stage IIIA"]<-3
#tmp[tmp="Stage IIIC"]<-3
#tmp[tmp="Stage IV"]<-4
#tmp[tmp="Stage IVA"]<-4
#tmp[tmp=="Stage I"]<-1
#curated$stageall<-tmp

tmp <- sapply(curated$unique_patient_ID,function(id)
  clinical.drug$drug_name[clinical.drug$bcr_patient_barcode==id])

tmp[sapply(tmp, length)==0] <- NA
tmp.therapy<-sapply(curated$unique_patient_ID,function(id)
    clinical.drug$therapy_type[clinical.drug$bcr_patient_barcode==id])
tmp.therapy[sapply(tmp.therapy, length)==0] <- NA

tmp.batch<-sapply(curated$unique_patient_ID,function(id)
	clinical.batch$BCR.batch[clinical.batch$bcr_patient_barcode==id])
tmp.batch[sapply(tmp.batch,length)==0]<- NA

curated$batch<-sapply(tmp.batch, function(x)
	return(x))
	

curated$fu <- sapply(tmp, function(x)
  ifelse(length(grep("5 FU | 5- FU | 5-Fluorouracil | 5-Fluoruouracil | 5FU | 5-FU | Fluorouracil", x))>0, "y",
         ifelse(is.na(x), NA,"n")))

curated$bevacizumab<- sapply(tmp, function(x) ifelse(length(grep("Bevacizumab | Avastin", x))>0,"y",
                                                     ifelse(is.na(x), NA, "n")))

curated$irinotecan<-sapply(tmp, function(x) ifelse(length(grep("Camptosar | Irinotecan | Irinotecan HCl", x))>0,"y",
                                                   ifelse(is.na(x), NA, "n")))

curated$capecitabine<-sapply(tmp, function(x) ifelse(length(grep("Capecitabine", x))>0,"y",
                                                     ifelse(is.na(x), NA, "n")))

curated$cpt11<-sapply(tmp, function(x) ifelse(length(grep("CPT-11", x))>0,"y",
                                              ifelse(is.na(x), NA, "n")))

curated$dexamethasone<-sapply(tmp, function(x) ifelse(length(grep("Dexamethasone", x))>0,"y",
                                                    ifelse(is.na(x), NA, "n")))
curated$erbitux<-sapply(tmp, function(x) ifelse(length(grep("Erbitux | Cetuximab", x))>0,"y",
                                                ifelse(is.na(x), NA, "n")))
curated$gcsf<-sapply(tmp, function(x) ifelse(length(grep("Filgrastim (G-CSF)", x))>0,"y",
                                             ifelse(is.na(x), NA, "n")))

curated$fudr<-sapply(tmp, function(x) ifelse(length(grep("Floxuridine", x))>0,"y",
                                             ifelse(is.na(x), NA, "n")))
curated$folfiri<-sapply(tmp, function(x) ifelse(length(grep("Folfiri", x))>0,"y",
                                                ifelse(is.na(x), NA, "n")))


curated$folfox<-sapply(tmp, function(x) ifelse(length(grep("Folfox | FOLFOX | Folfox-4", x))>0,"y",
                                               ifelse(is.na(x), NA, "n")))

curated$leucovorin<-sapply(tmp, function(x) ifelse(length(grep("Folinic acid | Leocovorin | Leucovorin | Levcovorin | Leucovorin Calcium", x))>0,"y",
                                            ifelse(is.na(x), NA, "n")))
curated$mitomycin<-sapply(tmp, function(x) ifelse(length(grep("Mitomycin | Mitomycin C", x))>0,"y",
                                                  ifelse(is.na(x), NA, "n")))
curated$platin<-sapply(tmp, function(x) ifelse(length(grep("Oxaliplatin | oxaliplatinum+ 5-FU", x))>0,"y",
                                               ifelse(is.na(x), NA, "n")))

curated$panitumumab<-sapply(tmp, function(x) ifelse(length(grep("Panitumumab", x))>0,"y",
                                                    ifelse(is.na(x), NA, "n")))

curated$pegfilgrastim<-sapply(tmp, function(x) ifelse(length(grep("Pegfilgrastim", x))>0, "y",
                                                      ifelse(is.na(x), NA, "n")))

curated$raltitrexed<-sapply(tmp, function(x) ifelse(length(grep("Raltitrexed", x))>0, "y",
                                                    ifelse(is.na(x), NA, "n")))
curated$xeloda<-sapply(tmp, function(x) ifelse(length(grep("XELODA | Xeloda", x))>0, "y",
                                               ifelse(is.na(x), NA, "n")))

curated$ancillary<-sapply(tmp.therapy, function(x) ifelse(length(grep("Ancillary", x))>0, "y",
                                                  ifelse(is.na(x), NA, "n")))

curated$chemotherapy<-sapply(tmp.therapy, function(x) ifelse(length(grep("Chemotherapy", x))>0, "y",
                                                  ifelse(is.na(x), NA, "n")))

curated$moltherapy<-sapply(tmp.therapy, function(x) ifelse(length(grep("Targeted Molecular therapy", x))>0, "y",
                                                  ifelse(is.na(x), NA, "n")))

curated <- postProcess(curated, uncurated, do.celfile.batch=FALSE)
