source setvars
mkdir ../FINAL
mkdir ../FINAL/DEFAULT
for file in `find $DATAHOME -name "*_curated_*"`;do rsync -v $file $REPOSHOME/FINAL/DEFAULT;done
for file in `find $DATAHOME -name "*rma_*"`;do rsync -v $file $REPOSHOME/FINAL/RMA;done
