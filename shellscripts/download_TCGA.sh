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

for file in `ls *.tar.gz`;do tar xfz $file;done
find . -name "*tcga_level3.data.txt" > TCGA_file_sources.txt
for file in `find . -name "*tcga_level3.data.txt"`;do mv $file .;done
rm -f *.tar.gz
