source setvars
mkdir $REPOSHOME/MAPPED

##Agilent

##$REXEC CMD BATCH --vanilla "--args GSE20842 $REPOSHOME/GENEMAPS/GSE20842.csv  $DATAHOME/GSE20842/PROCESSED/DEFAULT/GSE20842_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE20842_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE21815 $REPOSHOME/GENEMAPS/efg_agilent_wholegenome_4x44k_v1.csv $DATAHOME/GSE21815/PROCESSED/DEFAULT/GSE21815_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE21815_genemapping.log

##hg_u95av2

$REXEC CMD BATCH --vanilla "--args GSE11237 $REPOSHOME/GENEMAPS/affy_hg_u95av2.csv $DATAHOME/GSE11237/PROCESSED/RMA/GSE11237_rma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE11237_genemapping.log

##hg_u133_plus2

$REXEC CMD BATCH --vanilla "--args GSE13067 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE13067/PROCESSED/RMA/GSE13067_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE13067_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE13294 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE13294/PROCESSED/DEFAULT/GSE13294_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE13294_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE14095 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE14095/PROCESSED/DEFAULT/GSE14095_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE14095_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE14333 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE14333/PROCESSED/DEFAULT/GSE14333_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE14333_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE17536 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE17536/PROCESSED/DEFAULT/GSE17536_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE17536_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE17537 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE17537/PROCESSED/RMA/GSE17537_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE17537_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE17538-GPL570 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE17538-GPL570/PROCESSED/DEFAULT/GSE17538-GPL570_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE17538-GPL570_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE18105 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE18105/PROCESSED/RMA/GSE18105_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE18105_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE2109 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE2109/PROCESSED/DEFAULT/GSE2109_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE2109_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE21510 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE21510/PROCESSED/RMA/GSE21510_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE21510_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE26682-GPL570 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE26682-GPL570/PROCESSED/DEFAULT/GSE26682-GPL570_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE26682-GPL570_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE26906 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE26906/PROCESSED/RMA/GSE26906_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE26906_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE28702 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE28702/PROCESSED/RMA/GSE28702_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE28702_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE33113 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE33113/PROCESSED/RMA/GSE33113_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE33113_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE39582 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE39582/PROCESSED/DEFAULT/GSE39582_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE39582_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE4526 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE4526/PROCESSED/DEFAULT/GSE4526_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE4526_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE45270 $REPOSHOME/GENEMAPS/affy_hg_u133_plus_2.csv $DATAHOME/GSE45270/PROCESSED/RMA/GSE45270_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE45270_genemapping.log

##hg_u133a

$REXEC CMD BATCH --vanilla "--args GSE12945 $REPOSHOME/GENEMAPS/affy_hg_u133a.csv $DATAHOME/GSE12945/PROCESSED/RMA/GSE12945_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE12945_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE26682-GPL96 $REPOSHOME/GENEMAPS/affy_hg_u133a.csv $DATAHOME/GSE26682-GPL96/PROCESSED/RMA/GSE26682-GPL96_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE26682-GPL96_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE4045 $REPOSHOME/GENEMAPS/affy_hg_u133a.csv $DATAHOME/GSE4045/PROCESSED/RMA/GSE4045_frma_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE4045_genemapping.log

##huex-1_0-st

$REXEC CMD BATCH --vanilla "--args GSE16125-GPL5175 $REPOSHOME/GENEMAPS/affy_huex_1_0_st_v2.csv $DATAHOME/GSE16125-GPL5175/PROCESSED/DEFAULT/GSE16125-GPL5175_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE16125-GPL5175_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE24549-GPL5175 $REPOSHOME/GENEMAPS/affy_huex_1_0_st_v2.csv $DATAHOME/GSE24549-GPL5175/PROCESSED/DEFAULT/GSE24549-GPL5175_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE24549-GPL5175_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE24550-GPL5175 $REPOSHOME/GENEMAPS/affy_huex_1_0_st_v2.csv $DATAHOME/GSE24550-GPL5175/PROCESSED/DEFAULT/GSE24550-GPL5175_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE24550-GPL5175_genemapping.log

##custom
$REXEC CMD BATCH --vanilla "--args $REPOSHOME/GENEMAPS/GSE12225.csv $DATAHOME/GSE12225/PROCESSED/DEFAULT/GSE12225_default_gpl.csv GB_ACC embl" $SRCHOME/getPlatformsBiomaRt.R $LOG/GSE12225_platform.log
$REXEC CMD BATCH --vanilla "--args GSE12225 $REPOSHOME/GENEMAPS/GSE12225.csv $DATAHOME/GSE8842/PROCESSED/DEFAULT/GSE12225_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE12225_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE12520-GPL7199 $REPOSHOME/DATA/GSE12520-GPL7199/PROCESSED/DEFAULT/GSE12520-GPL7199_default_gpl.csv $DATAHOME/GSE12520-GPL7199/PROCESSED/DEFAULT/GSE12520-GPL7199_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE12520-GPL7199_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE12520-GPL9336 $REPOSHOME/DATA/GSE12520-GPL9336/PROCESSED/DEFAULT/GSE12520-GPL9336_default_gpl.csv $DATAHOME/GSE12520-GPL9336/PROCESSED/DEFAULT/GSE12520-GPL9336_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE12520-GPL9336_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE2630 $REPOSHOME/DATA/GSE2630/PROCESSED/DEFAULT/GSE2630_default_gpl.csv $DATAHOME/GSE2630/PROCESSED/DEFAULT/GSE2630_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE2630_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE27544 $REPOSHOME/DATA/GSE27544/PROCESSED/DEFAULT/GSE27544_default_gpl.csv $DATAHOME/GSE27544/PROCESSED/DEFAULT/GSE27544_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE27544_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE3294 $REPOSHOME/DATA/GSE3294/PROCESSED/DEFAULT/GSE3294_default_gpl.csv $DATAHOME/GSE3294/PROCESSED/DEFAULT/GSE3294_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE3294_genemapping.log
$REXEC CMD BATCH --vanilla "--args GSE3964 $REPOSHOME/DATA/GSE3964/PROCESSED/DEFAULT/GSE3964_default_gpl.csv $DATAHOME/GSE3964/PROCESSED/DEFAULT/GSE3964_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/GSE3964_genemapping.log

##TCGA
$REXEC CMD BATCH --vanilla "--args TCGA $REPOSHOME/GENEMAPS/TCGA.csv $DATAHOME/TCGA/PROCESSED/DEFAULT/TCGA_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/TCGA_genemapping.log
$REXEC CMD BATCH --vanilla "--args TCGA-RNASeqV2 $REPOSHOME/GENEMAPS/TCGA-RNASeqV2.csv $DATAHOME/TCGA-RNASeqV2/PROCESSED/DEFAULT/TCGA-RNASeqV2_default_curated_exprs.txt $REPOSHOME/MAPPED" $SRCHOME/geneMapper.R $LOG/TCGA-RNASeqV2_genemapping.log

