\name{ GSE2109_eset }
\alias{ GSE2109_eset }
\docType{data}
\title{ Expression Project for Oncology (expO) }
\description{ EXpression Project for Oncology, International Genomics Consortium, www.intgen.org }
\usage{data( GSE2109_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: EXpression Project for Oncology, International Genomics Consortium, www.intgen.org 
  Laboratory: expO, IGC 2005 
  Contact information:  
  Title: Expression Project for Oncology (expO) 
  URL:  
  PMIDs: PMID unknown 

  Abstract: A 8 word abstract is available. Use 'abstract' method.
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
assayData: 19320 features, 427 samples
Platform type: hgu133plus2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      427 character character 

sample_type: 
tumor 
  427 

primarysite: 
 co  re 
343  84 

summarygrade: 
high  low NA's 
  75  270   82 

G: 
   1    2    3    4 NA's 
  10  260   71    4   82 

summarystage: 
early  late  NA's 
  166   177    84 

T: 
   Length     Class      Mode 
      427 character character 

N: 
   0    1    2    X NA's 
 187  101   59    3   77 

M: 
   0    1    X NA's 
 274   64    9   80 

Dstage: 
   A    B    C    D NA's 
  51  113  105    5  153 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  25.00   55.00   65.00   65.68   75.00   95.00 

gender: 
   f    m NA's 
 207  219    1 

stageall: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.000   2.000   2.000   2.516   3.000   4.000      90 

ethnicity: 
    black caucasian  hispanic     other 
       23       391         1        12 

batch: 
   Length     Class      Mode 
      427 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
      427 character character 

}}
\keyword{datasets}
