platforms <- c("illumina_humanwg_6_v1",
               "illumina_humanwg_6_v2",
               "illumina_humanwg_6_v3",
               "illumina_humanht_12",
               "efg_agilent_wholegenome_4x44k_v1",
               "efg_agilent_wholegenome_4x44k_v2",
               "efg_agilent_sureprint_g3_ge_8x60k",
               "affy_hg_focus",
               "affy_hg_u133_plus_2",
               "affy_hg_u133a_2",
               "affy_hg_u133a",
               "affy_hg_u133b",
               "affy_hg_u95av2",
               "affy_hg_u95b",
               "affy_hg_u95c",
               "affy_hg_u95d",
               "affy_hg_u95e",
               "affy_hg_u95a",
               "affy_hugenefl",
	       "affy_huex_1_0_st_v2",
	       "affy_hugene_1_0_st_v1"
		)

library(biomaRt)
ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")

##listAttributes(ensembl)$name   #to see all attributes

for (platform in platforms){
  print(paste("Getting",platform))
  if(exists("genes")) rm(genes)
  genes <- getBM(attributes = c(platform, "hgnc_symbol", "ensembl_gene_id", "entrezgene"),mart=ensembl)
  genes <- genes[genes[,1]!="",]  #remove rows with no probeset ID
  genes <- genes[genes[,2]!="",]  #remove rows with no hgnc_symbol
  genes <- genes[,1:2]  #keep probeset and hgnc_symbol only
  genes <- unique(genes) #keep unique rows
  genemap <- tapply(genes[,2],genes[,1],paste,collapse="///")
  genemap <- data.frame(probeset=names(genemap),hgnc=genemap)
  fname <- paste(platform,".csv",sep="")
  write.csv(genemap,paste("../GENEMAPS/",fname,sep=""),row.names=FALSE,quote=FALSE)
}
