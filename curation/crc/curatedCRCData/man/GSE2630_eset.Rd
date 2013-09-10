\name{ GSE2630_eset }
\alias{ GSE2630_eset }
\docType{data}
\title{ A gene signature of 8 genes could identify the risk of recurrence and progression in Dukes' B colon cancer patients. }
\description{ The benefit of postoperative adjuvant chemotherapy in patients with Dukes' B colorectal cancer is still uncertain and its routine use is not recommended. The five-year relapse rate is approximately 25-40\% and the identification of patients at high risk of recurrence would represent an important strategy for the use of adjuvant chemotherapy. We retrospectively analyzed gene expression profiles in frozen tumor specimens from patients with Dukes' B colorectal cancer by using high density oligonucleotide microarrays. Our results show a subset of 48 genes differentially expressed with an associated probability <0.001 in the t-test. Another statistical procedure based on the Fisher criterion resulted in 11 genes able to separate both groups. We selected the 8 genes present in both subsets. The differential expression of five genes (CHD2, RPS5, ZNF148, BRI3 and MGC23401) in colon cancer progression was confirmed by real-time PCR in an independent set of patients of Dukes' B and C stages. }
\usage{data( GSE2630_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Bandr??E,??Malumbres R,??Cubedo E,??Sola J,??Garc??Foncillas J,??Labarga A. A gene signature of 8 genes could identify the risk of recurrence and progression in Dukes' B colon cancer patients. Oncol Rep 2007 May 
  Laboratory: Bandr??E,??Malumbres R,??Cubedo E,??Sola J,??Garc??Foncillas J,??Labarga A 2005 
  Contact information:  
  Title: A gene signature of 8 genes could identify the risk of recurrence and progression in Dukes' B colon cancer patients. 
  URL:  
  PMIDs: 17390049 

  Abstract: A 151 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Human 19K oligo array
   platform_shorttitle:     
      NA
   platform_summary:     
      human 19k oligo array
   platform_manufacturer:     
      NA
   platform_distribution:     
      non-commercial
   platform_accession:     
      GPL2006
   platform_technology:     
      spotted oligonucleotide
   warnings:     
      No warnings yet

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1 10 ... 9 (48 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 48 features, 16 samples
Platform type: human 19k oligo array
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       16 character character 

T: 
 2  3 
 6 10 

N: 
 0 
16 

M: 
 0 
16 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  42.00   59.50   67.00   64.38   72.00   81.00 

recurrence_status: 
norecurrence   recurrence 
          10            6 

summarylocation: 
   l    r NA's 
   2    4   10 

gender: 
 f  m 
 5 11 

uncurated_author_metadata: 
   Length     Class      Mode 
       16 character character 

}}
\keyword{datasets}
