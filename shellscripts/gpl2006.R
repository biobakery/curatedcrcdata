library(GEOquery)
platform <- getGEO('GPL2006',AnnotGPL=TRUE)
platform <-platform@dataTable
platform <- platform@table
colnames(platform) <- gsub(pattern=" ",replacement = '_',x=colnames(platform))
colnames(platform)
head(platform[,1:10])
genes<-platform$Gene_symbol
names(genes)<-platform$ID
genes<-genes[genes!=""]
genemap<-data.frame(probeset=names(genes),hgnc=genes)
fname<-"gpl-2006.csv"
write.csv(genemap,paste("../GENEMAPS/",fname,sep=""),row.names=FALSE,quote=FALSE)
