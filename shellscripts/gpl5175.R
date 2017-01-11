## GeneMap for GPL-5175 (Affy huex)


library(huex10sttranscriptcluster.db)
x <-  huex10sttranscriptclusterSYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
  # Get the SYMBOL for the first five probes
  xx[1:5]
  # Get the first one
  xx[[1]]
}

genes<-do.call(rbind.data.frame,xx)
genes <- unique(genes) #keep unique rows
genemap <- data.frame(probeset=rownames(genes),hgnc=genes)
colnames(genemap)<-c("probeset","hgnc")
fname <- "gpl-5175.csv"
write.csv(genemap,paste("../GENEMAPS/",fname,sep=""),row.names=FALSE,quote=FALSE)


