\name{ GSE16125.GPL5175_eset }
\alias{ GSE16125.GPL5175_eset }
\docType{data}
\title{ Integrative approach for prioritizing cancer genes in sporadic colon cancer. }
\description{ The current multistep carcinogenesis models of colon cancer do not fully capture the genetic heterogeneity of the disease, which is additionally complicated by the presence of passenger and driver genetic alterations. The aim of this study was to select in the context of this significant heterogeneity additional genes functionally related to colon cancer development. High-throughput copy number and gene expression data of 36 microsatellite stable sporadic colon cancers resected from patients of a single institution characterized for mutations in APC, KRAS, TP53 and loss of 18q were analyzed. Genes whose expression correlated with the underlying copy number pattern were selected, and their association with the above listed mutations and overall survival was evaluated. Gain of 20q was strongly associated with TP53 mutation, and overall survival with alterations on 7p, 8p, 13q, 18q, and 20q. An association with 18q loss and gain of 8q24 was also observed. New candidate genes with a potential role in colon cancer are PLCG1 on 20q, DBC1 on 8q21, and NDGR1 on 8p24. In addition, an unexpected pattern of loss and mutability was found in the region upstream of the KRAS gene. By integrating copy number alterations with gene expression and mutations in colon cancer associated genes, we have developed a strategy that identifies previously known molecular features and additional players in the molecular landscape of colon cancer.Copyright 2009 Wiley-Liss, Inc. }
\usage{data( GSE16125.GPL5175_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Reid JF, Gariboldi M, Sokolova V, Capobianco P et al.??Integrative approach for prioritizing cancer genes in sporadic colon cancer.??Genes Chromosomes Cancer??2009 Nov 
  Laboratory: Reid, Pierotti 2009 
  Contact information:  
  Title: Integrative approach for prioritizing cancer genes in sporadic colon cancer. 
  URL:  
  PMIDs: 19672874 

  Abstract: A 225 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array [transcript (gene) versio
n]
   platform_shorttitle:     
      Affymetrix HuEx-1_0-st
   platform_summary:     
      huex-1_0-st
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL5175
   platform_technology:     
      in situ oligonucleotide
   warnings:     
      No warnings yet

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: ABCF1 ANAPC5 ... ZNF259P1///ZNF259 (282 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 282 features, 36 samples
Platform type: huex-1_0-st
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   4 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 32.000  32.000  32.000  22.000   1.397   0.822      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       36 character character 

summarystage: 
early  late  NA's 
    6    27     3 

T: 
   1    2    3    4 NA's 
   1    5    8   19    3 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  19.00   61.75   67.50   65.58   71.25   90.00 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   30.0   202.5   390.0   926.2  1658.0  2970.0       4 

vital_status: 
deceased   living     NA's 
      22       10        4 

gender: 
 f  m 
20 16 

kras: 
 n  y 
16 20 

mutation_apc: 
 n  y 
19 17 

batch: 
   Length     Class      Mode 
       36 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
       36 character character 

}}
\keyword{datasets}
