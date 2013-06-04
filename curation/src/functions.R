mapstage <- function(x){
  output <- x
  output[output=="I"] <- 1
  output[output=="II"] <- 2
  output[output=="III"] <- 3
  output[output=="IV"] <- 4
  output <- as.integer(output)
  return(output)
}

initialCuratedDF <- function(DF.rownames,template.filename){
  template <- read.csv(template.filename,as.is=TRUE)
  output <- matrix(NA,
                   ncol=nrow(template),
                   nrow=length(DF.rownames))
  colnames(output) <- template$col.name
  rownames(output) <- DF.rownames
  output <- data.frame(output)
  for (i in 1:ncol(output)){
    class(output[,i]) <- template[i,"var.class"]
  }
  output$sample_name <- DF.rownames
  return(output)
}
