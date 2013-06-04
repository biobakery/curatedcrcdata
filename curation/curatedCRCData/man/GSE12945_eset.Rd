\name{GSE12945_eset}
\alias{GSE12945_eset}
\docType{data}
\title{An expression module of WIPF1-coexpressed genes identifies patients with favorable prognosis in three tumor types.}
\description{Wiskott-Aldrich syndrome (WAS) predisposes patients to leukemia and lymphoma. WAS is caused by mutations in the protein WASP which impair its interaction with the WIPF1 protein. Here, we aim to identify a module of WIPF1-coexpressed genes and to assess its use as a prognostic signature for colorectal cancer, glioma, and breast cancer patients. Two public colorectal cancer microarray data sets were used for discovery and validation of the WIPF1 co-expression module. Based on expression of the WIPF1 signature, we classified more than 400 additional tumors with microarray data from our own experiments or from publicly available data sets according to their WIPF1 signature expression. This allowed us to separate patient populations for colorectal cancers, breast cancers, and gliomas for which clinical characteristics like survival times and times to relapse were analyzed. Groups of colorectal cancer, breast cancer, and glioma patients with low expression of the WIPF1 co-expression module generally had a favorable prognosis. In addition, the majority of WIPF1 signature genes are individually correlated with disease outcome in different studies. Literature gene network analysis revealed that among WIPF1 co-expressed genes known direct transcriptional targets of c-myc, ESR1 and p53 are enriched. The mean expression profile of WIPF1 signature genes is correlated with the profile of a proliferation signature. The WIPF1 signature is the first microarray-based prognostic expression signature primarily developed for colorectal cancer that is instrumental in other tumor types: low expression of the WIPF1 module is associated with better prognosis.
We used microarrays for the validation of a WIPF1 co-expression module which was developed on two publically available datasets.}
\usage{data(GSE12945_eset)}
\format{

  Authors: Staub E, Groene J, Heinze M, Mennerich D et al.

  Lab: Signature Diagnostics AG

  Contact: bernd.hinzmann@signature-dignostics.de

  Title: An expression module of WIPF1-coexpressed genes identifies patients with favorable prognosis in three tumor types.

  Url: http://pubget.com/search?q=19399471 

  PubMedID: 19399471 
}
\details{
  Platform used: [HG-U133A] Affymetrix Human Genome U133A Array

  frma normalization used from the frma bioconductor package.

  62 samples.

  GEO_platform_accession: GSE12945
  platform_summary: GPL96
  biomart_ID: [HG-U133A] Affymetrix Human Genome U133A Array

  Overall survival time-to-event summary (in years):
  Call: survfit(formula = Surv(time, cens) ~ -1)

  ___________________________________________________________________

  Available sample meta-data


  ___________________________________________________________________

  alt_sample_name:
  Length = 62.
  Class = character.
  Mode = character.

  sample_type:
  tumor = 62.
  
  days_to_death: 
  Min. = 210 Median = 1395 Mean = 1266.77 Max. = 1920

}
\source{
  http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE12945
}
\keyword{datasets}