##01 - "Recurrent Solid Tumor" and 11 - "Solid Tissue Normal".
##histological_type -> sample_type
tmp <- tumor.num
tmp <- sub("01","tumor",tmp,fixed=TRUE)
tmp <- sub("11","adjacentnormal",tmp,fixed=TRUE)
curated$sample_type <- tmp

##primary_site -> tumor_tissue_site
tmp <- uncurated$tumor_tissue_site 
tmp[tmp=="Colon"] <- "co"
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
tmp[tmp=="Splenic Flexure"]<-NA
tmp[tmp=="Descending Colon"]<-"descending"
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

curated <- postProcess(curated, uncurated, do.celfile.batch=FALSE)