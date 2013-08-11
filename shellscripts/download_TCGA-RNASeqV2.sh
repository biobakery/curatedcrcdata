#!/bin/bash

source setvars

URL="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/coad/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/"
strInputAccession="TCGA-RNASeqV2"
strBaseDir=$DATAHOME
celDir=$strBaseDir/$strInputAccession
processedDir=$celDir/PROCESSED
defaultDir=$processedDir/DEFAULT

mkdir $strBaseDir
mkdir $celDir
mkdir $celDir/RAW
mkdir $processedDir
mkdir $defaultDir

cd $celDir/RAW

wget --timestamping --no-check-certificate $URL/unc.edu_COAD.IlluminaHiSeq_RNASeqV2.Level_3.1.4.0.tar.gz 
wget --timestamping --no-check-certificate $URL/unc.edu_COAD.IlluminaHiSeq_RNASeqV2.mage-tab.1.4.0.tar.gz  


for file in `ls *.tar.gz`;do tar xfz $file;done
find . -name "*.rsem.genes.normalized_results" > TCGA-RNASeqV2_file_sources.txt
for file in `find . -name "*.rsem.genes.normalized_results"`;do mv $file .;done
for file in `find . -name "*.sdrf.txt"`;do mv $file .;done
rm -f *.tar.gz

RCODE="
    library(reshape2);
    files = dir(full.names=TRUE);
    files = files[grep(\".rsem.genes.normalized_results\", files)];
     data <- lapply(files, read.delim, stringsAsFactors=FALSE, as.is=TRUE);
    sdrf <- read.delim(\"unc.edu_COAD.IlluminaHiSeq_RNASeqV2.1.3.0.sdrf.txt\",
        as.is=TRUE);
    barcodes <- sdrf[match(gsub(\"\\\\..*$\",\"\",gsub(\".*unc.edu.\",\"\",files)),
     sdrf[,1]), 2];
    data <- lapply(1:length(data), function(i) cbind(barcode=barcodes[i],
      data[[i]]));
    rdata <- do.call(rbind, data);
    cdata <- dcast(gene_id~barcode,data=rdata,
    value.var=\"normalized_count\");
    write.csv(cdata, file=\"TCGA-RNASeqV2_default_exprs.csv\",
    quote=FALSE, row.names=FALSE);
    hgnc <- gsub(\"\\\\|.*$\",\"\", cdata[,1]);
    write.csv(data.frame(probeset=cdata[,1],hgnc=hgnc),
    file=\"TCGA-RNASeqV2.csv\",
    quote=FALSE,row.names=FALSE);
    "

TMPFILE=`mktemp -t cod-rnaseqv2.XXXXXXXXXX`
echo $RCODE > "$TMPFILE"
R --vanilla < "$TMPFILE"
mv TCGA-RNASeqV2_default_exprs.csv $defaultDir
mkdir $REPOSHOME/GENEMAPS
mv TCGA-RNASeqV2.csv $REPOSHOME/GENEMAPS
mv unc.edu_COAD.IlluminaHiSeq_RNASeqV2.1.3.0.sdrf.txt $UNCURATED
