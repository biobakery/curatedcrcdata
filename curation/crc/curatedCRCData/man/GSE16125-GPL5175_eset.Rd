\name{GSE16125-GPL5175_eset}
\alias{GSE16125-GPL5175_eset}
\docType{data}
\title{Integrative approach for prioritizing cancer genes in sporadic colon cancer}
\description{The current multistep carcinogenesis models of colon cancer do not fully capture the genetic heterogeneity of the disease, which is additionally complicated by the presence of passenger and driver genetic alterations. The aim of the present study was to select in the context of this significant heterogeneity additional genes functionally related to colon cancer development. Methods: High-throughput copy number and gene expression data of 36 microsatellite stable sporadic colon cancers resected from patients of a single Institution characterized for mutations in APC, KRAS, TP53 and loss of 18q were analyzed. Genes whose expression correlated with the underlying copy number pattern were selected and their association with the above listed mutations and overall survival was evaluated. Results: Gain of 20q was strongly associated with TP53 mutation, and overall survival with alterations on 7p, 8p, 13q, 18q and 20q. An association with 18q loss and gain of 8q24 was also observed. New candidate genes with a potential role in colon cancer are PLCG1 on 20q, DBC1 on 8q21 and NDGR1 on 8p24. In addition an unexpected pattern of loss and mutability was found in the region upstream of KRAS gene. Conclusions: By integrating copy number alterations with gene expression and mutations in colon cancer associated genes we have developed a strategy that identifies previously known molecular features and additional players in the molecular landscape of colon cancer.}
\usage{data(GSE16125-GPL5175_eset)}
\format{

  Authors: Reid JF, Gariboldi M, Sokolova V, Capobianco P et al.

  Lab: Fondazione Istituto FIRC di Oncologia Molecolare

  Contact: james.reid@ifom-ieo-campus.it

  Title: Integrative approach for prioritizing cancer genes in sporadic colon cancer

  Url: http://pubget.com/search?q=19672874

  PubMedID: 19672874 
}
\details{
  Platform used: [HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]

  36 samples.

  GEO_platform_accession: GSE16125-GPL5175_eset
  platform_summary: GPL5175
  biomart_ID: [HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]

  Overall survival time-to-event summary (in years):
  Call: survfit(formula = Surv(time, cens) ~ -1)

  ___________________________________________________________________

  Available sample meta-data


  ___________________________________________________________________

  alt_sample_name:
  Length = 36.
  Class = character.
  Mode = character.

  sample_type:
  tumor = 36.
  
  days_to_death: 
  Min = 30 Median = 390 Mean = 926.25 Max = 2970

}
\source{
  http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE16125-GPL5175_eset
}
\keyword{datasets}