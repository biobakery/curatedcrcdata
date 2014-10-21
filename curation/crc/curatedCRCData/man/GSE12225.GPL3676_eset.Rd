\name{ GSE12225.GPL3676_eset }
\alias{ GSE12225.GPL3676_eset }
\docType{data}
\title{ Integrating chromosomal aberrations and gene expression profiles to dissect rectal tumorigenesis. }
\description{ Accurate staging of rectal tumors is essential for making the correct treatment choice. In a previous study, we found that loss of 17p, 18q and gain of 8q, 13q and 20q could distinguish adenoma from carcinoma tissue and that gain of 1q was related to lymph node metastasis. In order to find markers for tumor staging, we searched for candidate genes on these specific chromosomes.We performed gene expression microarray analysis on 79 rectal tumors and integrated these data with genomic data from the same sample series. We performed supervised analysis to find candidate genes on affected chromosomes and validated the results with qRT-PCR and immunohistochemistry.Integration of gene expression and chromosomal instability data revealed similarity between these two data types. Supervised analysis identified up-regulation of EFNA1 in cases with 1q gain, and EFNA1 expression was correlated with the expression of a target gene (VEGF). The BOP1 gene, involved in ribosome biogenesis and related to chromosomal instability, was over-expressed in cases with 8q gain. SMAD2 was the most down-regulated gene on 18q, and on 20q, STMN3 and TGIF2 were highly up-regulated. Immunohistochemistry for SMAD4 correlated with SMAD2 gene expression and 18q loss.On basis of integrative analysis this study identified one well known CRC gene (SMAD2) and several other genes (EFNA1, BOP1, TGIF2 and STMN3) that possibly could be used for rectal cancer characterization. }
\usage{data( GSE12225.GPL3676_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Lips EH, van Eijk R, de Graaf EJ, Oosting J et al.??Integrating chromosomal aberrations and gene expression profiles to dissect rectal tumorigenesis.BMC Cancer??2008 Oct 29 
  Laboratory: Lips,  Morreau 2008 
  Contact information:  
  Title: Integrating chromosomal aberrations and gene expression profiles to dissect rectal tumorigenesis. 
  URL:  
  PMIDs: 18959792 

  Abstract: A 221 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      NKI-CMF Homo sapiens 35k oligo array
   platform_shorttitle:     
      NA
   platform_summary:     
      nki-cmf homo sapiens 35k oligo array
   platform_manufacturer:     
      Central Microarray Facility, NKI Amsterdam
   platform_distribution:     
      non-commercial
   platform_accession:     
      GPL3676
   platform_technology:     
      spotted oligonucleotide
   warnings:     
      No warnings yet

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1CF ... ZZZ3 (19727 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19727 features, 42 samples
Platform type: nki-cmf homo sapiens 35k oligo array
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       42 character character 

sample_type: 
tumor 
   42 

primarysite: 
re 
42 

summarystage: 
early  late 
   30    12 

T: 
 1  2  3 
14 27  1 

N: 
 0  1  2 
30  6  6 

M: 
 0 
42 

location: 
rectum 
    42 

summarylocation: 
 l 
42 

stageall: 
 1  3 
29 13 

preop_drug_treatment: 
 n 
42 

uncurated_author_metadata: 
   Length     Class      Mode 
       42 character character 

}}
\keyword{datasets}
