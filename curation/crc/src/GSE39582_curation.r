
rm(list=ls())

source("../../functions.R")

uncurated1 <- read.csv("../uncurated/GSE39582_full_pdata.csv",as.is=TRUE,row.names=1)
uncurated2 <- read.csv("../uncurated/GSE_39582 addtl data.csv",as.is=TRUE)
uncurated <- merge( uncurated1, uncurated2, by.x="title", by.y="id")

celfile.dir <- "../../../DATA/GSE39582/RAW"
rownames(uncurated)<-uncurated$geo_accession
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_crc.csv")

##title -> alt_sample_name
curated$alt_sample_name <- uncurated$geo_accession


##sample_type
curated$sample_type<-"tumor"

##primarysite
curated$primarysite<-"co"

##N
tmp1<-uncurated$characteristics_ch1.3
tmp1[tmp1=="tnm.stage: 4"]<-NA
tmp1[tmp1=="tnm.stage: 2"]<-"0"
tmp1[tmp1=="tnm.stage: 1"]<-"0"
tmp1[tmp1=="tnm.stage: 3"]<-NA
tmp1[tmp1=="tnm.stage: 0"]<-"0"
curated$N<-tmp1

#M
tmp2<-uncurated$characteristics_ch1.3
tmp2[tmp2=="tnm.stage: 4"]<-"1"
tmp2[tmp2=="tnm.stage: 2"]<-"0"
tmp2[tmp2=="tnm.stage: 1"]<-"0"
tmp2[tmp2=="tnm.stage: 3"]<-"0"
tmp2[tmp2=="tnm.stage: 0"]<-"0"
curated$M<-tmp2

##gender
tmp<-apply(uncurated,1,getVal,string="Sex: ")
tmp<-sub("Sex: ","",tmp)
tmp[tmp=="M"]<-"m"
tmp[tmp=="F"]<-"f"
curated$gender<-tmp

##age_at_initial_pathologic_diagnosis
tmp<-uncurated$characteristics_ch1.2
tmp<-sub("age.at.diagnosis: ","",tmp)
tmp<-round(as.numeric(tmp))
curated$age_at_initial_pathologic_diagnosis<-tmp

##stageall
tmp<-apply(uncurated,1,getVal,string="tnm.stage: ")
tmp<-sub("tnm.stage: ","",tmp)
curated$stageall<-tmp

##summarystage
tmp<-uncurated$characteristics_ch1.3
tmp<-sub("tnm.stage: ","",tmp)
tmp[tmp=="2"]<-NA
tmp[tmp=="3"]<-"late"
tmp[tmp=="4"]<-"late"
tmp[tmp=="1"]<-"early"
tmp[tmp=="0"]<-"early"
curated$summarystage<-tmp

##Dstage
tmp3<-uncurated$characteristics_ch1.3
tmp3[tmp3=="tnm.stage: 4"]<-"D"
tmp3[tmp3=="tnm.stage: 2"]<-"B"
tmp3[tmp3=="tnm.stage: 1"]<-NA
tmp3[tmp3=="tnm.stage: 3"]<-"C"
tmp3[tmp3=="tnm.stage: 0"]<-NA
curated$Dstage<-tmp3

##location
tmp<-uncurated$characteristics_ch1.4
tmp<-sub("tumor.location: ","",tmp)
curated$location<-tmp

##summarylocation
tmp[tmp=="distal"]<-"l"
tmp[tmp=="proximal"]<-"r"
curated$summarylocation<-tmp

##kras
tmpkras<-uncurated$KRAS.Mutation  
tmpkras[tmpkras=="M"]<-"mutant"
tmpkras[tmpkras=="WT"]<-"wt"
curated$kras<-tmpkras  

##braf
tmpbraf<-uncurated$BRAF.Mutation
tmpbraf[tmpbraf=="M"]<-"mutant"
tmpbraf[tmpbraf=="WT"]<-"wt"
curated$braf<-tmpbraf  

##preop_drug_treatment
curated$preop_drug_treatment<-"n"

##dfs_status
tmp<-uncurated$characteristics_ch1.6
tmp<-sub("rfs.event: ","",tmp)
tmp[tmp==0]<-"living_norecurrence"
tmp[tmp==1]<-"deceased_or_recurrence"
curated$dfs_status<-tmp

##days_to_recurrence_or_death
tmp<-uncurated$characteristics_ch1.7  
tmp<-sub("rfs.delay: ","",tmp)
tmp <- as.numeric(tmp)
tmp <- tmp * 30  #months to days
curated$days_to_recurrence_or_death<-tmp

##chemotherapy
tmp<-uncurated$characteristics_ch1.5
tmp<-sub("chemotherapy.adjuvant: ","",tmp)
tmp[tmp=="N"]<-"n"
tmp[tmp=="Y"]<-"y"
curated$chemotherapy<-tmp

##msi
tmp<-uncurated$characteristics_ch1.8
tmp<-sub("mmr.status: ","",tmp)
tmp[tmp=="pMMR"]<-"MSS"
tmp[tmp=="dMMR"]<-"MSI"
curated$msi<-tmp

##drug_treatment
tmp<-uncurated$characteristics_ch1.5
tmp<-sub("chemotherapy.adjuvant: ","",tmp)
tmp[tmp=="N"]<-"n"
tmp[tmp=="Y"]<-"y"
curated$drug_treatment<-tmp

curated <- postProcess(curated, uncurated)

curated<-updatedfs(curated)

write.table(curated, row.names=FALSE, file="../curated/GSE39582_curated_pdata.txt",sep="\t")

