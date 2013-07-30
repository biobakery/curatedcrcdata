source setvars

cd $CURATIONSRC

rm $CURATED/*.txt
rm $CURATED/ERRORS/*

for file in `ls *.r`;do $REXEC CMD BATCH --vanilla $file $LOG/$file.log;done
$REXEC CMD BATCH --vanilla "--args $CURATED $CURATIONSRC/template_crc.csv" \
    $SRCHOME/checkCurated.R $LOG/checkCurated.log

cp $LOG/checkCurated.log $CURATED/ERRORS/

cd $CURATIONSRC/misc
$REXEC CMD BATCH --vanilla summarizeCurated.R $LOG/summarizeCurated.log
