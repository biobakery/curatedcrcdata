#!/bin/bash

source setvars

URL="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/read/cgcc/unc.edu/agilentg4502a_07_3/transcriptome/"

strInputAccession="TCGA-READ"
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

wget --timestamping $URL/unc.edu_READ.AgilentG4502A_07_3.Level_3.1.1.0.tar.gz
wget --timestamping $URL/unc.edu_READ.AgilentG4502A_07_3.Level_3.2.0.0.tar.gz
wget --timestamping $URL/unc.edu_READ.AgilentG4502A_07_3.mage-tab.1.3.0.tar.gz


for file in `ls *.tar.gz`;do tar xfz $file;done
find . -name "*tcga_level3.data.txt" > TCGA_file_sources.txt
for file in `find . -name "*tcga_level3.data.txt"`;do mv $file .;done
for file in `find . -name "*.sdrf.txt"`;do mv $file .;done
rm -f *.tar.gz

sdrf="unc.edu_READ.AgilentG4502A_07_3.sdrf.txt"

$REXEC CMD BATCH --vanilla "--args READ $sdrf" $REPOSHOME/shellscripts/tcga_PROCESSING.R $REPOSHOME/LOG/tcga_PROCESSING.log

mv TCGA-READ_default_exprs.csv $defaultDir
mkdir $REPOSHOME/GENEMAPS
mv TCGA-READ.csv $REPOSHOME/GENEMAPS
mv unc.edu_READ.AgilentG4502A_07_3.sdrf.txt $UNCURATED

