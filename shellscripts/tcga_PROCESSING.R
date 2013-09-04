#inputargs<-c("sample","sdrf")
inputargs<- commandArgs(TRUE)
sample<-inputargs[1]
file.sdrf<-inputargs[2]
library(reshape2)
files = dir(full.names=TRUE)
files = files[grep("tcga_level3.data.txt", files)]
    data <- lapply(files, read.delim, stringsAsFactors=FALSE, as.is=TRUE,skip=1)
                   sdrf <- read.delim(paste("unc.edu_",sample,".AgilentG4502A_07_3.sdrf.txt",sep=""),
        as.is=TRUE)
                   sdrf<-sdrf[which(sdrf[,1]!="Stratagene Univeral Reference"),]
    barcodes <- sdrf[match(gsub("\\\\..*$","",gsub("./","",files)),
     sdrf[,40]), 15]
                   data <- lapply(1:length(data), function(i) cbind(barcode=barcodes[i],
                   data[[i]]))
                   rdata <- do.call(rbind, data)
                   cdata <- dcast(Composite.Element.REF~barcode,data=rdata,
                   value.var="log2.lowess.normalized..cy5.cy3..collapsed.by.gene.symbol")
    write.csv(cdata, file=paste("TCGA-",sample,"_default_exprs.csv",sep=""),
    quote=FALSE, row.names=FALSE)
                   hgnc <- cdata[,1]
                   write.csv(data.frame(probeset=cdata[,1],hgnc=hgnc),
                   file=paste("TCGA-",sample,".csv",sep=""),
    quote=FALSE,row.names=FALSE)
