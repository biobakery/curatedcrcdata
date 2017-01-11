\name{ GSE26906_eset }
\alias{ GSE26906_eset }
\docType{data}
\title{ Expression Profiles in Stage II Colon Cancer According to APC Gene Status. }
\description{ Colorectal cancer is one of the most common cancers in the world. Histoclinical staging is efficient, but combination with molecular markers may improve the classification of stage II cancers. Several tumor-suppressor genes have been associated with colorectal cancer, and the most frequent allelic losses have been extensively studied for their prognosis effect, but the results remain controversial. In a previous study, we found a possible influence of the chromosome 5 status in the development of liver metastases in stage II colon cancers. We have here investigated the role of the APC gene, located in chromosome arm 5q, in a series of 183 colon adenocarcinomas through a combined analysis of gene expression, mutation, allelic loss and promoter methylation, and metastasis occurrence. Point mutations were found in 73\% of cases and allelic losses were found in 39\%; 59\% of tumors presented with a biallelic inactivation, with a very strong interdependence of the two APC hits (P = 2.1 x 10(-9)). No association was found between expression, number and type of APC alterations, and metastatic evolution. Our results show that the determination of APC status cannot help in the prediction of metastasis and cannot be used to subclassify stage II colon cancers. }
\usage{data( GSE26906_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Birnbaum DJ, Laibe S, Ferrari A, Lagarde A et al.??Expression Profiles in Stage II Colon Cancer According to APC Gene Status.??Transl Oncol??2012 Apr 
  Laboratory: Olschwang S 2011 
  Contact information:  
  Title: Expression Profiles in Stage II Colon Cancer According to APC Gene Status. 
  URL:  
  PMIDs: 22496922 

  Abstract: A 199 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
   platform_shorttitle:     
      Affymetrix HG-U133Plus2
   platform_summary:     
      hgu133plus2
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL570
   platform_technology:     
      in situ oligonucleotide
   warnings:     
      No warnings yet

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1BG-AS1 ... ZZZ3 (19320 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19320 features, 90 samples
Platform type: hgu133plus2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       90 character character 

sample_type: 
tumor 
   90 

primarysite: 
co 
90 

N: 
 0 
90 

M: 
 0 
90 

Dstage: 
 B 
90 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  31.00   58.50   69.50   66.52   75.00   94.00 

recurrence_status: 
norecurrence   recurrence 
          69           21 

msi: 
MSS 
 90 

summarylocation: 
 l  r 
65 25 

gender: 
 f  m 
47 43 

mutation_apc: 
 n  y 
22 68 

stageall: 
 2 
90 

batch: 
   Length     Class      Mode 
       90 character character 

drug_treatment: 
 n 
90 

preop_drug_treatment: 
 n 
90 

chemotherapy: 
 n 
90 

uncurated_author_metadata: 
   Length     Class      Mode 
       90 character character 

}}
\keyword{datasets}
