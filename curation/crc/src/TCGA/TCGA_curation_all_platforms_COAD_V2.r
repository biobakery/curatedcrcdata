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

##ethnicity
race<-uncurated$race
race[race=="WHITE"]<-"caucasian"
race[race=="BLACK OR AFRICAN AMERICAN"]<-"black"
curated$ethnicity<-race

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

##recurrence_status
recur<-uncurated$new_tumor_event_after_initial_treatment
recur[recur=="NO"]<-"norecurrence"
recur[recur=="YES"]<-"recurrence"
curated$recurrence_status<- recur

##days_to_tumor_recurrence  
curated$days_to_tumor_recurrence<-uncurated$days_to_new_tumor_event_after_initial_treatment

##MSI
tmp<-uncurated$microsatellite_instability
tmp[tmp=="YES"]<-"MSI"
tmp[tmp=="NO"]<-"MSS"
curated$msi<-tmp

##location
tmp<-uncurated$anatomic_neoplasm_subdivision
tmp[tmp=="Sigmoid Colon"]<-"sigmoid"
tmp[tmp=="Transverse Colon"]<-"transverse"
tmp[tmp=="Cecum"]<-"caecum"
tmp[tmp=="Ascending Colon"]<-"ascending"
tmp[tmp=="Hepatic Flexure"]<-"hepaticflexure"
tmp[tmp=="Splenic Flexure"]<-"splenicflexure"
tmp[tmp=="Descending Colon"]<-"descending"
tmp[tmp=="Rectosigmoid Junction"]<-"rectosigmoid" 
tmp[tmp=="Rectum"]<-"rectum"
curated$location<-tmp

##summarylocation
tmp<-uncurated$anatomic_neoplasm_subdivision
tmp[tmp=="Sigmoid Colon"]<-"l"
tmp[tmp=="Transverse Colon"]<-"r"
tmp[tmp=="Cecum"]<-"r"
tmp[tmp=="Ascending Colon"]<-"r"
tmp[tmp=="Hepatic Flexure"]<-"r"
tmp[tmp=="Splenic Flexure"]<-"l"
tmp[tmp=="Descending Colon"]<-"l"
tmp[tmp=="Rectosigmoid Junction"]<-"l" 
tmp[tmp=="Rectum"]<-"l"
curated$summarylocation<-tmp

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

##T
tmp<-uncurated$pathologic_T
tmp[tmp=="Tis"]<-"is"
tmp[tmp=="T0"]<-0
tmp[tmp=="T1"]<-1
tmp[tmp=="T2"]<-2
tmp[tmp=="T3"]<-3
tmp[tmp=="T4"]<-4
tmp[tmp=="T4a"]<-4
tmp[tmp=="T4b"]<-4
curated$T<-tmp

##N
tmp2<-uncurated$pathologic_N
tmp2[tmp2=="N0"]<-0
tmp2[tmp2=="N1"]<-1
tmp2[tmp2=="N1a"]<-1
tmp2[tmp2=="N1b"]<-1
tmp2[tmp2=="N1c"]<-1
tmp2[tmp2=="N2"]<-2
tmp2[tmp2=="N2a"]<-2
tmp2[tmp2=="N2b"]<-2
tmp2[tmp2=="NX"]<-"X"
curated$N<-tmp2

##M
tmp<-uncurated$pathologic_M
tmp[tmp=="M0"]<-0
tmp[tmp=="M1a"]<-1
tmp[tmp=="M1b"]<-1
tmp[tmp=="M1"]<-1
tmp[tmp=="MX"]<-"X"
curated$M<-tmp

# ##Filling in the version of AJCC classification used
# tmp<-uncurated$system_version
# tmp[1]<-"6th"
# tmp[12]<-"6th"
# tmp[29]<-"6th"
# tmp[43]<-"6th"
# tmp[53]<-"6th"
# tmp[109]<-"6th"
# tmp[110]<-"6th"
# tmp[124]<-"6th"
# tmp[127]<-"6th"
# tmp[128]<-"6th"
# uncurated$system_version<-tmp

##summarystage
summarystage<-tmp
tmp1 <-curated$T
tmp2 <-curated$N
tmp3 <-curated$M
summarystage[(tmp1<4) & (tmp2==0) & (tmp3==0)]<-"early"
summarystage[(tmp2>0) | (tmp3>0)] <-"late"
summarystage[(tmp1==4) & (tmp2==0) & (tmp3==0)]<-"late"
curated$summarystage <-summarystage

##stageall
tmp<-uncurated$pathologic_stage
tmp[tmp=="Stage IA"]<-1  
tmp[tmp=="Stage IIA"]<-2
tmp[tmp=="Stage II"]<-2
tmp[tmp=="Stage IIB"]<-2
tmp[tmp=="Stage IIIB"]<-3
tmp[tmp=="Stage III"]<-3
tmp[tmp=="Stage IIIA"]<-3
tmp[tmp=="Stage IIIC"]<-3
tmp[tmp=="Stage IV"]<-4
tmp[tmp=="Stage IVA"]<-4
tmp[tmp=="Stage I"]<-1
tmp[tmp=="NA"]<-NA
curated$stageall<-tmp

##lymphnodesremoved  
curated$lymphnodesremoved <- uncurated$lymph_node_examined_count

##lymphnodesinvaded
curated$lymphnodesinvaded<- uncurated$number_of_lymphnodes_positive_by_he

##preop_drug_treatment
neoadj<-uncurated$history_of_neoadjuvant_treatment
neoadj[neoadj=="No"]<-"n"
neoadj[neoadj=="Yes"]<-"y"
curated$preop_drug_treatment<-neoadj

##drug_treatment
drugrx<-uncurated$additional_pharmaceutical_therapy  
drugrx[drugrx=="NO"]<-"n"
drugrx[drugrx=="YES"]<-"y"
curated$drug_treatment<-drugrx

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
  ifelse(length(grep("5 FU | 5- FU | 5-Fluorouracil | 5-Fluoruouracil | 5FU | 5-FU | Fluorouracil | fluorouracil | 
                     oxaliplatinum+ 5-FU | Calcium Foliatum, fluorouracilum, oxaliplatinum, dexamethassone", x))>0, "y",
         ifelse(is.na(x), NA,"n")))

curated$bevacizumab<- sapply(tmp, function(x) ifelse(length(grep("Bevacizumab | Avastin", x))>0,"y",
                                                     ifelse(is.na(x), NA, "n")))

curated$irinotecan<-sapply(tmp, function(x) ifelse(length(grep("Camptosar | Irinotecan | Irinotecan HCl |IRINOtecan HCl | 
                                                               CPT-11 | Irinotecan HCL", x))>0,"y",
                                                   ifelse(is.na(x), NA, "n")))

curated$capecitabine<-sapply(tmp, function(x) ifelse(length(grep("Capecitabine | XELODA | Xeloda", x))>0,"y",
                                                     ifelse(is.na(x), NA, "n")))

curated$dexamethasone<-sapply(tmp, function(x) ifelse(length(grep("Dexamethasone | Calcium Foliatum, fluorouracilum, oxaliplatinum, dexamethassone", x))>0,"y",
                                                    ifelse(is.na(x), NA, "n")))

curated$cetuximab<-sapply(tmp, function(x) ifelse(length(grep("Erbitux | Cetuximab | Cetuximab Study drug", x))>0,"y",
                                                ifelse(is.na(x), NA, "n")))

curated$gcsf<-sapply(tmp, function(x) ifelse(length(grep("Filgrastim (G-CSF)", x))>0,"y",
                                             ifelse(is.na(x), NA, "n")))

curated$fudr<-sapply(tmp, function(x) ifelse(length(grep("Floxuridine", x))>0,"y",
                                             ifelse(is.na(x), NA, "n")))

curated$folfiri<-sapply(tmp, function(x) ifelse(length(grep("Folfiri", x))>0,"y",
                                                ifelse(is.na(x), NA, "n")))


curated$folfox<-sapply(tmp, function(x) ifelse(length(grep("Folfox | FOLFOX | Folfox-4 | FolFox", x))>0,"y",
                                               ifelse(is.na(x), NA, "n")))

curated$leucovorin<-sapply(tmp, function(x) ifelse(length(grep("Folinic acid | Leocovorin | Leucovorin | Levcovorin | Leucovorin Calcium", x))>0,"y",
                                            ifelse(is.na(x), NA, "n")))

curated$mitomycin<-sapply(tmp, function(x) ifelse(length(grep("Mitomycin | Mitomycin C", x))>0,"y",
                                                  ifelse(is.na(x), NA, "n")))

curated$platin<-sapply(tmp, function(x) ifelse(length(grep("Oxaliplatin | oxaliplatinum+ 5-FU | Calcium Foliatum, fluorouracilum, 
                                                           oxaliplatinum, dexamethassone | oxaliplatin", x))>0,"y",
                                               ifelse(is.na(x), NA, "n")))

curated$panitumumab<-sapply(tmp, function(x) ifelse(length(grep("Panitumumab", x))>0,"y",
                                                    ifelse(is.na(x), NA, "n")))

curated$pegfilgrastim<-sapply(tmp, function(x) ifelse(length(grep("Pegfilgrastim | Pegfilgrastim (Peg G-CSF))", x))>0, "y",
                                                      ifelse(is.na(x), NA, "n")))

curated$raltitrexed<-sapply(tmp, function(x) ifelse(length(grep("Raltitrexed", x))>0, "y",
                                                    ifelse(is.na(x), NA, "n")))

curated$ancillary<-sapply(tmp.therapy, function(x) ifelse(length(grep("Ancillary", x))>0, "y",
                                                  ifelse(is.na(x), NA, "n")))

curated$chemotherapy<-sapply(tmp.therapy, function(x) ifelse(length(grep("Chemotherapy", x))>0, "y",
                                                  ifelse(is.na(x), NA, "n")))

curated$moltherapy<-sapply(tmp.therapy, function(x) ifelse(length(grep("Targeted Molecular therapy", x))>0, "y",
                                                  ifelse(is.na(x), NA, "n")))

curated <- postProcess(curated, uncurated, do.celfile.batch=FALSE)
