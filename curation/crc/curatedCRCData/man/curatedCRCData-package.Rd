\name{curatedCRCData-package}
\alias{curatedCRCData}
\docType{package}
\title{Clinically Annotated Data for the CRC Cancer Transcriptome}
\description{The curatedCRCData package provides manually curated
  clinical data, uniformly processed expression data, and convenience
  functions for gene expression analysis in patients with colorectal cancer.
  
  Note: For staging, the "summarystage" was curated as using information on the pathologic T, N, and M stages, such that "early" = T stage of Tis through T3 and N0 and M0 and "late" = T stage of 4 and/or any N > 0 or M > 0}
\details{
\tabular{ll}{
Package: \tab curatedCRCData\cr
Type: \tab Package\cr
Version: \tab \Sexpr{packageDescription("curatedCRCData")$Version}\cr
Date: \tab \Sexpr{packageDescription("curatedCRCData")$Date}\cr
License: \tab Artistic-2.0\cr
Depends: \tab R (>= 2.10.0), affy\cr
}
}
\author{
Princy Parsana, Markus Riester, Levi Waldron

Department of Biostatistics and Computational Biology, Dana-Farber Cancer Institute, Harvard School of Public Health

Maintainer: Princy Parsana <princyparsana@jhu.edu>
}

\examples{
##List all datasets:
data(package="curatedCRCData")
##
}
