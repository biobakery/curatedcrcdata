\name{ GSE4045_eset }
\alias{ GSE4045_eset }
\docType{data}
\title{ Serrated carcinomas form a subclass of colorectal cancer with distinct molecular basis. }
\description{ Serrated colorectal carcinomas (CRCs) are morphologically different from conventional CRCs and have been proposed to follow a distinct pathway of CRC formation. Despite studies of single molecular events in this tumor type, the diagnosis of serrated CRC relies on morphology and the putative unique biological character of these tumors has not been established. Here we show that the gene expression profiling of 37 CRCs separated serrated and conventional CRCs into two distinct branches in unsupervised hierarchical clustering (P-value 7.8 x 10(-7)), and revealed 201 differentially expressed genes representing potential biomarkers for serrated CRC. Immunohistochemistry was utilized to verify the key findings in the 37 CRCs examined by expression profiling, and a separate validation set of 37 serrated and 86 conventional CRCs was examined to evaluate the candidate biomarkers in an extended sample material. Ephrin receptor B2, hypoxia-inducible factor 1-alpha and patched appeared as proteins important for genesis of serrated CRC. This study establishes serrated CRCs as a biologically distinct subclass of CRC and represents a step forward in the molecular classification of these cancers. The study also provides a platform to understand the molecular basis of serrated CRC and in long term may contribute to the development of specific treatment options for this tumor type. }
\usage{data( GSE4045_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Laiho P, Kokko A, Vanharanta S, Salovaara R et al.??Serrated carcinomas form a subclass of colorectal cancer with distinct molecular basis.??Oncogene2007 Jan 11 
  Laboratory: Laiho, Aaltonen 2007 
  Contact information:  
  Title: Serrated carcinomas form a subclass of colorectal cancer with distinct molecular basis. 
  URL:  
  PMIDs: 16819509 

  Abstract: A 205 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG-U133A] Affymetrix Human Genome U133A Array
   platform_shorttitle:     
      Affymetrix HG-U133A
   platform_summary:     
      hgu133a
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL96
   platform_technology:     
      in situ oligonucleotide
   warnings:     
      No warnings yet

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (12986 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 12986 features, 37 samples
Platform type: hgu133a
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       37 character character 

sample_type: 
tumor 
   37 

summarygrade: 
 low NA's 
  36    1 

G: 
   1    2    3 NA's 
   4   28    4    1 

summarystage: 
early  late 
    2    35 

M: 
 0  1 
35  2 

Dstage: 
 B  C  D 
 2 33  2 

family_history: 
   n    y NA's 
  32    3    2 

msi: 
MSI MSS 
  1  36 

location: 
  distal proximal 
      22       15 

summarylocation: 
 l  r 
22 15 

gender: 
 f  m 
19 18 

stageall: 
   3    4 NA's 
  33    2    2 

batch: 
   Length     Class      Mode 
       37 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
       37 character character 

}}
\keyword{datasets}
