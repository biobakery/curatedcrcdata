original <- read.delim("../uncurated/clinical_patient_coad.txt",as.is=TRUE,row.names=1,na.strings=c("[Not Available]", "[Not Applicable]", "[Pending]"),head=T)

##follow-up data:
follow.up <- read.delim("../uncurated/clinical_follow_up_v1.0_coad.txt",as.is=TRUE,row.names=1,na.strings=c("[Not Available]", "[Not Applicable]", "[Pending]"))

##stromal data:
clinical.slide <-
  read.delim("../uncurated/biospecimen_aliquot_coad.txt", as.is=TRUE,
             na.strings=c("[Not Available]", "[Not Reported]", "[Not Applicable]", "[Pending]"))
clinical.slide$bcr_sample_barcode <-
  as.factor(clinical.slide$bcr_sample_barcode)

##treatment data:
clinical.drug <-
  read.delim("../uncurated/clinical_drug_coad.txt", as.is=TRUE,
             na.strings=c("[Not Available]", "[Not Reported]", "[Not Applicable]", "[Pending]"))

clinical.drug$bcr_patient_barcode <-
  as.factor(clinical.drug$bcr_patient_barcode)

patient.IDs <- sapply( strsplit(rownames(follow.up), split="-"), function(x) paste(x[1:3], collapse="-") )
##Split the follow.up dataframe into a list, one element per unique patient ID:
follow.up.list <- lapply( unique(patient.IDs), function(id)
  follow.up[grep(id, rownames(follow.up)), ] )
##order each list element by year, month, day of form completion, then
##keep the first element.  The first element will either be the only
##element, or the one occuring at the *latest* date:
follow.up.list <- lapply(follow.up.list, function(x){
  x[order(x$date_of_form_completion,
          decreasing=TRUE)[1], ]
})
##Now rbind back to a single dataframe, which contains only the most recent follow-up:
follow.up <- do.call(rbind, follow.up.list)

##Get rid of what comes after the third "-" to change rownames to patient IDs:
rownames(follow.up) <- sapply( strsplit(rownames(follow.up), split="-"), function(x) paste(x[1:3], collapse="-") )

##Now we need to merge the original and the follow-up dataframes.
##First, get rid of columns in the original dataframe which are also
##present in the follow-up dataframe:
original <- original[, !colnames(original) %in% colnames(follow.up)]

##Some patients were present in the original, but not in
##the follow.up.  Add rows of NA to these patients in the follow-up:
missing.in.followup <- rownames(original)[!rownames(original) %in% rownames(follow.up)]
blank.df <- data.frame(matrix(NA, nrow=length(missing.in.followup), ncol=ncol(follow.up)))
rownames(blank.df) <- missing.in.followup
colnames(blank.df) <- colnames(follow.up)
##Then merge these blank rows with the original:
follow.up <- rbind(follow.up, blank.df)

##Finally merge the original and follow-up clinical data:
if( all(rownames(follow.up) %in% rownames(original)) & all(rownames(original) %in% rownames(follow.up)) )
  follow.up <- follow.up[match(rownames(original), rownames(follow.up)), ]
if( identical(rownames(original), rownames(follow.up)) )
  uncurated <- cbind(original, follow.up)


