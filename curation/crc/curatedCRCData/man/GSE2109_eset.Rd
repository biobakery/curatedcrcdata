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
      hg-u133_plus_2
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
assayData: 19320 features, 428 samples
Platform type: hg-u133_plus_2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      428 character character 

summarygrade: 
high  low NA's 
  75  270   83 

G: 
   1    2    3    4 NA's 
  10  260   71    4   83 

summarystage: 
early  late  NA's 
  166   168    94 

T: 
   1    2    3    4 NA's 
  14   55  244   31   84 

N: 
   0    1    2 NA's 
 187  101   59   81 

M: 
   0    1 NA's 
 274   64   90 

family_history: 
   n    y NA's 
 184  241    3 

gender: 
   f    m NA's 
 208  219    1 

stageall: 
   1    2    3    4 NA's 
  52  114  104   64   94 

ethnicity: 
other 
  428 

uncurated_author_metadata: 
   Length     Class      Mode 
      428 character character 

}}
\keyword{datasets}
