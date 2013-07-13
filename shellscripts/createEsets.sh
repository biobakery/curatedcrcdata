source setvars
$REXEC CMD BATCH --vanilla "--args ../MAPPED/DATA MaxMean ../etc/crc_cancer_public_datasets.csv" $SRCHOME/createEsets.R $LOG/createEsets_maxMean.log
$REXEC CMD BATCH --vanilla "--args ../MAPPED/DATA Normalizer ../etc/crc_cancer_public_datasets.csv" $SRCHOME/createEsets.R $LOG/createEsets_Normalizer.log
$REXEC CMD BATCH --vanilla "--args ../MAPPED/DATA FULL ../etc/crc_cancer_public_datasets.csv" $SRCHOME/createEsets.R $LOG/createEsets_FULL.log
