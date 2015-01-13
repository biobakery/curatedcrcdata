##This script should reproduce the entire repository, assuming curation is already completed and NCBI blast is installed (See README.md in main directory).

source setvars
mkdir $REPOSHOME/LOG
mkdir $REPOSHOME/DATA

#Install needed packages:
$REXEC CMD BATCH --vanilla $SRCHOME/install_needed_packages.R $LOG/install_needed_packages.log

./downloadRAW.sh
./downloadPROCESSED.sh
./curation.sh       #only to re-run the curation scripts and re-test the output.
./attach_PROCESSED.sh
./preprocess_AFFY.sh

##Create gene maps from Blast and BioMART.  May be run simultaneously:
$REXEC CMD BATCH --vanilla $SRCHOME/blast_gene_maps.R $LOG/blast_gene_maps.log
$REXEC CMD BATCH --vanilla $SRCHOME/getPlatforms.R $LOG/getPlatforms.log


#Do the mapping.
$REXEC CMD BATCH --vanilla gpl5175.R $LOG/gpl5175.log
$REXEC CMD BATCH --vanilla gpl2006.R $LOG/gpl2006.log
$REXEC CMD BATCH --vanilla gpl2829.R $LOG/gpl2829.log
./crcgenemapper.sh
$REXEC CMD BATCH --vanilla $SRCHOME/collapseRows.R $LOG/collapseRows.log

#Move the unmapped and two mapped versions to their own directories.
mkdir $REPOSHOME/MAPPED/DATA
mkdir $REPOSHOME/MAPPED/DATA/FULL/
mkdir $REPOSHOME/MAPPED/DATA/MaxMean/
mkdir $REPOSHOME/MAPPED/DATA/MaxMean_probesused/
mkdir $REPOSHOME/MAPPED/DATA/Normalizer/
mv $REPOSHOME/MAPPED/*FULL.txt $REPOSHOME/MAPPED/DATA/FULL/
mv $REPOSHOME/MAPPED/*probesused_MaxMean.txt $REPOSHOME/MAPPED/DATA/MaxMean_probesused/
mv $REPOSHOME/MAPPED/*MaxMean.txt $REPOSHOME/MAPPED/DATA/MaxMean/
mv $REPOSHOME/MAPPED/*Normalizer.txt $REPOSHOME/MAPPED/DATA/Normalizer/

#Copy the metadata into one place:
mkdir $REPOSHOME/MAPPED/DATA/metadata/
rm $REPOSHOME/MAPPED/DATA/metadata/*
find $DATAHOME -name "*_curated_pdata.txt" | xargs cp --target-dir=$REPOSHOME/MAPPED/DATA/metadata/
# that works on Macs
find $DATAHOME -name "*_curated_pdata.txt" | xargs -I % cp % $REPOSHOME/MAPPED/DATA/metadata/

#Finally, create esets:
./createEsets.sh


##Create help pages:
$REXEC CMD BATCH "--args $REPOSHOME/MAPPED/DATA MaxMean" $SRCHOME/summarizeEset.R $LOG/summarizeEset.Rout
$REXEC CMD BATCH "--args $REPOSHOME/MAPPED/DATA FULL" $SRCHOME/summarizeEset.R $LOG/summarizeEsetFULL.Rout
$REXEC CMD BATCH "--args $REPOSHOME/MAPPED/DATA Normalizer" $SRCHOME/summarizeEset.R $LOG/summarizeEsetNormalizer.Rout


#create curatedcrcdata package:
./curatedcrcdata.sh
