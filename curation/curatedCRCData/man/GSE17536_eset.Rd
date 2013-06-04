\name{GSE17536_eset}
\alias{GSE17536_eset}
\docType{data}
\title{Experimentally derived metastasis gene expression profile predicts recurrence and death in patients with colon cancer.}
\description{Staging inadequately predicts metastatic risk in colon cancer patients. We used a gene expression profile derived from invasive murine colon cancer cells that were highly metastatic in an immunocompetent mouse model to identify colon cancer patients at risk for recurrence in a phase I, exploratory biomarker study.}
\usage{data(GSE17536_eset)}
\format{

  Authors: Smith JJ, Deane NG, Wu F, Merchant NB et al.

  Lab: Vanderbilt University

  Contact: daniel.beauchamp@vanderbilt.edu

  Title: Integrative approach for prioritizing cancer genes in sporadic colon cancer

  Url: http://pubget.com/search?q=19914252 

  PubMedID: 19914252  
}
\details{
  Platform used: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

  frma normalization used.

  177 samples.

  GEO_platform_accession: GSE17536_eset
  platform_summary: GPL570
  biomart_ID: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

  Overall survival time-to-event summary (in years):
  Call: survfit(formula = Surv(time, cens) ~ -1)

  ___________________________________________________________________

  Available sample meta-data


  ___________________________________________________________________

  alt_sample_name:
  Length = 177.
  Class = character.
  Mode = character.

  sample_type:
  tumor = 177.
  
  days_to_death: 
  Min = 27.6 Median = 1268.1 Mean = 1444.67 Max = 4276.5

}
\source{
  http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE17536_eset
}
\keyword{datasets}