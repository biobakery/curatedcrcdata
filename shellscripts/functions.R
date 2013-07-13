##Function for merging replicate samples using a specified function
FUN_mergeDups <-
  function(obj,sourcename,FUN="mean")
{
  ##obj is a matrix, ExpressionSet, or LumiBatch object
  ##
  ##if obj is a p x n matrix, sourcename is a vector giving the source
  ##of each column.  Columns with the same source will be grouped
  ##together.  If obj is an ExpressionSet of LumiBatch object,
  ##sourcename is a name of one of the columns of pdata(obj), which
  ##will be used as the sample source for merging.
  ##
  ##FUN is the character name of a function to be applied to
  ##duplicates, such that get(FUN) returns the desired function.  If
  ##FUN="mean", FUN_mergeDups uses rowMeans to speed up computation.
  if(class(sourcename) != "character" & class(sourcename) != "factor" & class(sourcename) != "integer") stop("sourcename should be character, factor, or integer")
  if(class(obj) == "ExpressionSet")
    {
      library(affy)
      vec <- pData(obj)[[sourcename]]  #replicate info
      dat <- exprs(obj)                #matrix of numeric data
    }else if(class(obj) == "LumiBatch"){
      library(lumi)
      vec <- pData(obj)[[sourcename]]
      dat <- exprs(obj)
    }else if(class(obj) == "matrix" | class(obj) == "dataframe"){
      dat <- as.matrix(obj)
      vec <- sourcename
      if(length(vec) != nrow(dat)) stop("If obj is a matrix or dataframe, vec should have the same length as nrow(obj)")
    }else{
      stop("class(obj) should be ExpressionSet, LumiBatch, matrix, or dataframe")
    }
  dat.list <- lapply(unique(vec),function(x) dat[,vec==x])
  names(dat.list) <- make.names(unique(vec))
  myFUN <- get(FUN)
  dat.merged <-
    sapply(dat.list,function(x)
           {
             if(class(x)=="matrix"){
               if(FUN=="mean"){
                 ##do a trick for speed-up
                 return(rowMeans(x))
               }else{
                 return(apply(x,1,myFUN))
               }
             }else{
               return(x)
             }
           }
           )
  if(class(obj) == "ExpressionSet" | class(obj) == "LumiBatch"){
    pdata <- pData(obj)
    pdata.list <- lapply(unique(vec),function(x) pdata[vec==x,])

    pdata.merged <- sapply(pdata.list,function(pdata.one){
        apply(pdata.one,2,function(x)
          {
              if( all(is.na(x)) )
                  return(NA)
              if( length(unique(x)) == 1)
                  return( unique(x) )
###If the field separator "///" is in x, then merge within each field.
###This should only happen for uncurated_author_metadata.
              if(length(grep("///", x)) > 1){
                  x.split <- strsplit(x, split="///")
                  if(identical(length(x.split[[1]]), length(x.split[[2]]))){
###Loop over each element defined by "///" (this 
                      output <- sapply(1:length(x.split[[1]]), function(i){
                          x.split.unique <- unique(c(x.split[[1]][i], x.split[[2]][i]))
                          if(length(x.split.unique) == 1)
                              return(x.split.unique)
                          x.split.unique.matrix <- do.call(rbind, strsplit(x.split.unique, split=": "))
                          if(length(unique(x.split.unique.matrix[,1])) == 1){
                              return(paste(x.split.unique.matrix[1, 1], ": ", paste(x.split.unique.matrix[, 2], collapse=" / "), sep="") )
                          }else{
                              return(paste(x.split.unique, collapse=" / ") )
                          }
                      }
                                       )
                      output <- paste(output, collapse="///")
                  }
                  return( output )
              }else{
###Otherwise, if metadata for technical replicates conflict return NA:
                  return( NA )
              }
              ## return(paste(unique(na.omit(x)),collapse="/"))
          }
              )
    }
                           )

    pdata.merged <- t(pdata.merged)
    pdata.merged <- data.frame(pdata.merged)
    rownames(pdata.merged) <- make.names(unique(vec))
    pdata.merged <- pdata.merged[match(colnames(dat.merged),rownames(pdata.merged)),]
    pdata.merged.adf <- AnnotatedDataFrame(data=pdata.merged)
    if(identical(all.equal(sampleNames(pdata.merged.adf),colnames(dat.merged)),TRUE))
      {
        obj.merged <- new(class(obj),
                          phenoData=pdata.merged.adf,
                          exprs=dat.merged,
                          annotation=obj@annotation)
        dat.merged <- obj.merged
      }
  }
  return(dat.merged)
}
         
