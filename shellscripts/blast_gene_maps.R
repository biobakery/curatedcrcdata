##blast executable:
##ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.25+-x64-linux.tar.gz

blast.url <- "ftp://ftp.ensembl.org/pub/release-68/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.68.cdna.all.fa.gz"
blast.target <- basename(blast.url)

genemap.dir <- "../GENEMAPS"
dir.create(genemap.dir)

##Get a map between Ensembl gene ID and HGNC gene symbol
library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensg.hgnc <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = ensembl)
##ensg.hgnc <- ensg.hgnc[ensg.hgnc$hgnc_symbol!="",]

setwd(genemap.dir)

##get the target cDNA fasta file from Ensembl:
if(!file.exists(paste(genemap.dir,blast.target,sep="/"))){
  download.file(blast.url, destfile=blast.target)
}

system(paste("gzip -cd ",blast.target," > tmp.fa"))  #un-gzip
system("sed -r -i -e 's/ cdna.+gene:/_/g' tmp.fa")   #Change names to ENSTnnnnnn_ENSGnnnnn
system(paste("makeblastdb -in tmp.fa -out ", sub(".fa.gz", "", blast.target, fixed=TRUE), " -dbtype nucl", sep=""))  #make the database
unlink("tmp.fa")  #clean up



##Make the target file when sequences are available in the default_gpl.csv, so we can later blast these sequences
for (studyid in dir("../DATA")){
  print(studyid)
  platform.file <- paste("../DATA/",studyid,"/PROCESSED/DEFAULT/",studyid,"_default_gpl.csv",sep="")
  if(!file.exists(platform.file)) next
  platform <- read.csv(platform.file,as.is=TRUE)
  oligo.count <- sapply(platform[1:1000,],function(x) sum(grepl("^[ACGT]{25,}$",x)))  ##Count how many of each column look like an oligo sequence, ie 25 or more ACGT, and nothing else.
  if(any(oligo.count > 500)){  #only require half because sometimes the oligos are missing for some probes
    sequence.col <- which.max(oligo.count)
    make.fasta <- TRUE
    print(paste("The column titled:",names(sequence.col)))
    print("looks like oligo sequences.")
  }else{
    make.fasta <- FALSE
    print("Nothing in the platform file looks like oligo sequences.")
  }
  if(make.fasta){
    filename <- paste(genemap.dir,"/",studyid,"_target.fa",sep="")
    unlink(filename)
    for (i in 1:nrow(platform)){
      this.sequence <- platform[i,sequence.col]
      if(grepl("^[ACGT]{25,}$",this.sequence)){  #if this looks like an oligo sequence, add it to the FASTA file:
        write.table(paste(">",platform[i,1],sep=""),file=filename,row.names=FALSE,col.names=FALSE,append=TRUE,quote=FALSE)
        write.table(this.sequence,file=filename,row.names=FALSE,col.names=FALSE,append=TRUE,quote=FALSE)
      }
    }
  }
}


##Do the blast:
for (this.target in dir(pattern="_target.fa")){
  this.id <- sub("_target.fa","",this.target)
  ##do the blast
  system(paste("blastn -query ", this.target,
               " -db ", sub(".fa.gz", "", blast.target, fixed=TRUE),
               " -outfmt 6 -evalue 1e-10 -word_size 20 > ",this.id,"_blastresults.txt",sep=""))
  #Read the blast results
  genes <- read.delim(paste(this.id,"_blastresults.txt",sep=""),as.is=TRUE,header=FALSE)
  genes <- genes[,1:2]  #just keep probeset ID and ENST_ENSG ID
  genes[,2] <- sub("^ENST[0-9]+_","",genes[,2])  #get rid of transcript IDs, just keep gene IDs
  genes$hgnc <- ensg.hgnc[match(genes[,2],ensg.hgnc[,1]),2]  ##add hgnc IDs
  genes <- genes[,c(1,3)]  #just keep probe ID and hgnc
  genes <- unique(genes)
  genes <- genes[genes$hgnc!="",]
  genemap <- tapply(genes[,2],genes[,1],paste,collapse="///")  #genes[,2] is the symbol
  genemap <- data.frame(probeset=names(genemap),hgnc=genemap)
  platform <- sub("_target.fa","",this.target,fixed=TRUE)
  fname <- paste(platform,".csv",sep="")
  write.csv(genemap,paste("../GENEMAPS/",fname,sep=""),row.names=FALSE,quote=FALSE)
}














