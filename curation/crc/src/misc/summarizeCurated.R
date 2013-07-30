CURATED <- "../../curated"
cwd <- getwd()
setwd(CURATED)
allfiles <- dir(pattern="*.txt")

allstudies <- lapply(allfiles,function(thisfile){
  x <- read.delim(thisfile,as.is=TRUE)
  x$study <- sub("_curated_pdata.txt","",thisfile)
  return(x)
})

names(allstudies) <- allfiles

sort(sapply(allstudies,ncol))

allstudies <- do.call(rbind,allstudies)


write.csv(allstudies,"crc.studies.summary.csv")