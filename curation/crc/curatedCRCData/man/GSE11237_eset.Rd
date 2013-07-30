\name{GSE11237_eset}
\alias{GSE11237_eset}
\docType{data}
\title{Celecoxib pre-treatment in human colorectal adenocarcinoma
  patients.}
\description{A cohort of 23 patients with primary colorectal
  adenocarcinomas was randomized to receive a 7-day course of celecoxib
  (400 mg b.i.d) or no drug prior to surgical resection.  Gene
  expression profiling was performed on resected adenocarcinomas from
  patients with and without celecoxib pre-treatment.  Using fold change
  (>1.5) and p-value (<0.05) cut-offs, 190 genes were differentially
  expressed between adenocarcinomas from patients receiving celecoxib
  and those that did not...Celecoxib pre-treatment for 7 days in vivo is
  associated with alterations in colorectal adenocarcinoma gene
  expression which are suggestive of diminished cellular proliferation}
\usage{data(GSE11237_eset)}
\format{

  Authors: Auman JT, Church R, Lee SY, Watson MA et al.

  Lab: Institute for Pharmacogenomics and Individualized Therapy at the
  University of North Carolina at Chapel Hill.

  Contact: jtauman@email.unc.edu

  Title: Celecoxib pre-treatment in human colorectal adenocarcinoma patients.

  Url: http://pubget.com/search?q=18653328

  PubMedID: 18653328
}
\details{
  Platform used: [HG_U95Av2] Affymetrix Human Genome U95 Version 2
  Array.

  rma normalization used from the rma bioconductor package.

  23 samples.

  GEO_platform_accession: GSE11237
  platform_summary: GPL8300
  biomart_ID: HG_U95Av2 affy_hg_u95av2

  Overall survival time-to-event summary (in years):
  Call: survfit(formula = Surv(time, cens) ~ -1)

  ___________________________________________________________________

  Available sample meta-data


  ___________________________________________________________________

  alt_sample_name:
  Length = 23.
  Class = character.
  Mode = character.

  sample_type:
  tumor = 23.

}
\source{
  http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE11237
}
\keyword{datasets}