#!/bin/bash

source setvars

URL="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/coad/cgcc/unc.edu/agilentg4502a_07_3/transcriptome/"

strInputAccession="TCGA"
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

wget --timestamping $URL/unc.edu_COAD.AgilentG4502A_07_3.Level_3.1.5.0.tar.gz 
wget --timestamping $URL/unc.edu_COAD.AgilentG4502A_07_3.Level_3.2.0.0.tar.gz
wget --timestamping $URL/unc.edu_COAD.AgilentG4502A_07_3.mage-tab.1.9.0.tar.gz


for file in `ls *.tar.gz`;do tar xfz $file;done
find . -name "*tcga_level3.data.txt" > TCGA_file_sources.txt
for file in `find . -name "*tcga_level3.data.txt"`;do mv $file .;done
for file in `find . -name "*.sdrf.txt"`;do mv $file .;done
rm -f *.tar.gz

RCODE="
    library(reshape2);
    files = dir(full.names=TRUE);	
    files = files[grep(\"tcga_level3.data.txt\", files)];
    data <- lapply(files, read.delim, stringsAsFactors=FALSE, as.is=TRUE,skip=1);
    sdrf <- read.delim(\"unc.edu_COAD.AgilentG4502A_07_3.sdrf.txt\",
        as.is=TRUE);
    sdrf<-sdrf[which(sdrf[,\"Source.Name\"]!=\"Stratagene Univeral Reference\"),];
    barcodes <- sdrf[match(gsub(\"\\\\..*$\",\"\",gsub(\"./\",\"\",files)),
     sdrf[,40]), 15];
    data <- lapply(1:length(data), function(i) cbind(barcode=barcodes[i],
      data[[i]]));
    rdata <- do.call(rbind, data);
     cdata <- dcast(Composite.Element.REF~barcode,data=rdata,
     value.var=\"log2.lowess.normalized..cy5.cy3..collapsed.by.gene.symbol\");
    write.csv(cdata, file=\"TCGA_default_exprs.csv\",
    quote=FALSE, row.names=FALSE);
    hgnc <- cdata[,1];
    write.csv(data.frame(probeset=cdata[,1],hgnc=hgnc),
    file=\"TCGA.csv\",
    quote=FALSE,row.names=FALSE);
    "

TMPFILE=`mktemp -t cod-tcga.XXXXXXXXXX`
echo $RCODE > "$TMPFILE"
R --vanilla < "$TMPFILE"
mv TCGA_default_exprs.csv $defaultDir
mkdir $REPOSHOME/GENEMAPS
mv TCGA.csv $REPOSHOME/GENEMAPS
mv unc.edu_COAD.AgilentG4502A_07_3.sdrf.txt $UNCURATED
