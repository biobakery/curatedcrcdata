##This script copies ExpressionSets into a shell R package directory
##and builds the MaxMean, Normalizer, and FULL versions of the
##curatedOvarianData R package.

source setvars

cd $REPOSHOME/curation/crc/
rm *.tar.gz
rm *.pdf
rm -r $REPOSHOME/curation/crc/curatedCRCData/data
rm -r NormalizerVcuratedCRCData
rm -r FULLVcuratedCRCData
mkdir $REPOSHOME/curation/crc/curatedCRCData/data

##prepare directories for Normalizer version
cp -a $REPOSHOME/curation/crc/curatedCRCData $REPOSHOME/curation/crc/NormalizerVcuratedCRCData
sed -i -e 's/curatedCRCData/NormalizerVcuratedCRCData/g' $REPOSHOME/curation/crc/NormalizerVcuratedCRCData/DESCRIPTION
sed -i -e 's/curatedCRCData/NormalizerVcuratedCRCData/g' $REPOSHOME/curation/crc/NormalizerVcuratedCRCData/inst/extdata/*
mv $REPOSHOME/curation/crc/NormalizerVcuratedCRCData/man/curatedCRCData-package.Rd $REPOSHOME/curation/crc/NormalizerVcuratedCRCData/man/NormalizerVcuratedCRCData-package.Rd
sed -i -e 's/curatedCRCData/NormalizerVcuratedCRCData/g' $REPOSHOME/curation/crc/NormalizerVcuratedCRCData/man/NormalizerVcuratedCRCData-package.Rd
sed -i -e 's/curatedCRCData/NormalizerVcuratedCRCData/g' $REPOSHOME/curation/crc/NormalizerVcuratedCRCData/inst/unitTests/*
sed -i -e 's/curatedCRCData/NormalizerVcuratedCRCData/g' $REPOSHOME/curation/crc/NormalizerVcuratedCRCData/tests/*
cp $REPOSHOME/MAPPED/DATA/Normalizer/*.rda $REPOSHOME/curation/crc/NormalizerVcuratedCRCData/data
cp $REPOSHOME/MAPPED/DATA/Normalizer/*.Rd $REPOSHOME/curation/crc/NormalizerVcuratedCRCData/man
rm -r $REPOSHOME/curation/crc/NormalizerVcuratedCRCData/inst/doc

##prepare directories for FULL version:
cp -a $REPOSHOME/curation/crc/curatedCRCData $REPOSHOME/curation/crc/FULLVcuratedCRCData
sed -i -e 's/curatedCRCData/FULLVcuratedCRCData/g' $REPOSHOME/curation/crc/FULLVcuratedCRCData/DESCRIPTION
sed -i -e 's/curatedCRCData/FULLVcuratedCRCData/g' $REPOSHOME/curation/ovarian/FULLVcuratedOvarianData/inst/extdata/*
mv $REPOSHOME/curation/crc/FULLVcuratedCRCData/man/curatedCRCData-package.Rd $REPOSHOME/curation/crc/FULLVcuratedCRCData/man/FULLVcuratedCRCData-package.Rd
sed -i -e 's/curatedCRCData/FULLVcuratedCRCData/g' $REPOSHOME/curation/crc/FULLVcuratedCRCData/man/FULLVcuratedCRCData-package.Rd
sed -i -e 's/curatedCRCData/FULLVcuratedCRCData/g' $REPOSHOME/curation/crc/FULLVcuratedCRCData/inst/unitTests/*
sed -i -e 's/curatedCRCData/FULLVcuratedCRCData/g' $REPOSHOME/curation/crc/FULLVcuratedCRCData/tests/*
cp $REPOSHOME/MAPPED/DATA/FULL/*.rda $REPOSHOME/curation/crc/FULLVcuratedCRCData/data
cp $REPOSHOME/MAPPED/DATA/FULL/*.Rd $REPOSHOME/curation/crc/FULLVcuratedCRCData/man
rm -r $REPOSHOME/curation/crc/FULLVcuratedCRCData/inst/doc


##copy data to MaxMean:
cp $REPOSHOME/MAPPED/DATA/MaxMean/*.rda $REPOSHOME/curation/crc/curatedCRCData/data
cp $REPOSHOME/MAPPED/DATA/MaxMean/*.Rd $REPOSHOME/curation/crc/curatedCRCData/man

##build the packages
$REXEC CMD build curatedCRCData
$REXEC CMD build NormalizerVcuratedCRCData
$REXEC CMD build FULLVcuratedCRCData
$REXEC CMD check *.tar.gz
$REXEC CMD Rd2pdf curatedCRCData
$REXEC CMD Rd2pdf NormalizerVcuratedCRCData
$REXEC CMD Rd2pdf FULLVcuratedCRCData

##clean-up
# rm $REPOSHOME/curation/crc/curatedCRCData/data/*.rda
# rm -rf $REPOSHOME/curation/crc/NormalizerVcuratedCRCData
# rm -rf $REPOSHOME/curation/crc/FULLVcuratedCRCData
