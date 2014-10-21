\name{ GSE11237_eset }
\alias{ GSE11237_eset }
\docType{data}
\title{ Celecoxib pre-treatment in human colorectal adenocarcinoma patients is associated with gene expression alterations suggestive of diminished cellular proliferation. }
\description{ Cancer cells treated with the cyclooxygenase-2 inhibitor celecoxib show growth inhibition and induced apoptosis. This study was conducted to determine if the same processes are relevant to celecoxib's effects on human colorectal adenocarcinomas treated in vivo. A cohort of 23 patients with primary colorectal adenocarcinomas was randomised to receive a 7-d course of celecoxib (400mg b.i.d.) or no drug prior to surgical resection. Gene expression profiling was performed on resected adenocarcinomas from the cohort of patients. Using fold change (>1.5) and p-value (<0.05) cut-offs, 190 genes were differentially expressed between adenocarcinomas from patients receiving celecoxib and those that did not. The celecoxib pre-treated samples showed decreased expression levels in multiple genes involved in cellular lipid and glutathione metabolism; changes associated with diminished cellular proliferation. Celecoxib pre-treatment for 7 d in vivo is associated with alterations in colorectal adenocarcinoma gene expression which are suggestive of diminished cellular proliferation. }
\usage{data( GSE11237_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Auman JT, Church R, Lee SY, Watson MA, Fleshman JW, Mcleod HL.??Celecoxib pre-treatment in human colorectal adenocarcinoma patients is associated with gene expression alterations suggestive of diminished cellular proliferation.??Eur J Cancer??2008 Aug 
  Laboratory: Auman, Mcleod 2008 
  Contact information:  
  Title: Celecoxib pre-treatment in human colorectal adenocarcinoma patients is associated with gene expression alterations suggestive of diminished cellular proliferation. 
  URL:  
  PMIDs: 18653328 

  Abstract: A 147 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Affymetrix Human Genome U95 Version 2 Array
   platform_shorttitle:     
      Affymetrix HG_U95Av2
   platform_summary:     
      hgu95av2
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL8300
   platform_technology:     
      in situ oligonucleotide
   warnings:     
      No warnings yet

Preprocessing: rma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: AADAC AAK1 ... ZZZ3 (8933 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 8933 features, 23 samples
Platform type: hgu95av2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       23 character character 

sample_type: 
tumor 
   23 

summarygrade: 
high  low 
   4   19 

G: 
 1  2  3 
 3 16  4 

summarystage: 
early  late 
   12    11 

T: 
 1  2  3  4 
 1  3 16  3 

N: 
 0  1  2 
15  7  1 

M: 
 0  1 
17  6 

location: 
   Length     Class      Mode 
       23 character character 

summarylocation: 
 l  r 
 9 14 

gender: 
 f  m 
13 10 

stageall: 
1 2 3 4 
4 9 4 6 

batch: 
2003-08-19 2003-08-28 
        15          8 

preop_drug_treatment: 
 n  y 
12 11 

preop_drug_name: 
celecoxib      NA's 
       11        12 

uncurated_author_metadata: 
   Length     Class      Mode 
       23 character character 

}}
\keyword{datasets}
