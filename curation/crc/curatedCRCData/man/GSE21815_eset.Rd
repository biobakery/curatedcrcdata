\name{ GSE21815_eset }
\alias{ GSE21815_eset }
\docType{data}
\title{ Long noncoding RNA HOTAIR regulates polycomb-dependent chromatin modification and is associated with poor prognosis in colorectal cancers. }
\description{ The functional impact of recently discovered long noncoding RNAs (ncRNAs) in human cancer remains to be clarified. One long ncRNA which has attracted attention is the Hox transcript antisense intergenic RNA termed HOTAIR, a long ncRNA expressed from the developmental HOXC locus located on chromosome 12q13.13. In cooperation with Polycomb complex PRC2, the HOTAIR long ncRNA is reported to reprogram chromatin organization and promote breast cancer metastasis. In this study, we examined the status and function of HOTAIR in patients with stage IV colorectal cancer (CRC) who have liver metastases and a poor prognosis. HOTAIR expression levels were higher in cancerous tissues than in corresponding noncancerous tissues and high HOTAIR expression correlated tightly with the presence of liver metastasis. Moreover, patients with high HOTAIR expression had a relatively poorer prognosis. In a subset of 32 CRC specimens, gene set enrichment analysis using cDNA array data revealed a close correlation between expression of HOTAIR and members of the PRC2 complex (SUZ12, EZH2, and H3K27me3). Our findings suggest that HOTAIR expression is associated with a genome-wide reprogramming of PRC2 function not only in breast cancer but also in CRC, where upregulation of this long ncRNA may be a critical element in metastatic progression. }
\usage{data( GSE21815_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Kogo R, Shimamura T, Mimori K, Kawahara K et al.??Long noncoding RNA HOTAIR regulates polycomb-dependent chromatin modification and is associated with poor prognosis in colorectal cancers.??Cancer Res??2011 Oct 15 
  Laboratory: Mori M,??Mimori K,??Yokobori T 2010 
  Contact information:  
  Title: Long noncoding RNA HOTAIR regulates polycomb-dependent chromatin modification and is associated with poor prognosis in colorectal cancers. 
  URL:  
  PMIDs: 21862635 

  Abstract: A 201 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name vers
ion)
   platform_shorttitle:     
      Agilent G4112F
   platform_summary:     
      agilent-014850 whole human genome microarray 4x44k g4112f
   platform_manufacturer:     
      Agilent
   platform_distribution:     
      commercial
   platform_accession:     
      GPL6480
   platform_technology:     
      in situ oligonucleotide
   warnings:     
      No warnings yet

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1BG-AS1 ... ZZZ3 (19686 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19686 features, 141 samples
Platform type: agilent-014850 whole human genome microarray 4x44k g4112f
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      141 character character 

sample_type: 
adjacentnormal          tumor           NA's 
             5             66             70 

summarystage: 
early  late  NA's 
   14    52    75 

T: 
   1    2    3    4 NA's 
   4   10   31   21   75 

N: 
   0    1    2 NA's 
  43   15    8   75 

M: 
   0    1 NA's 
  56   10   75 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  32.00   57.00   65.00   64.73   74.00   85.00       4 

gender: 
   f    m NA's 
  59   78    4 

stageall: 
   1    2    3    4 NA's 
  13   26   17   10   75 

uncurated_author_metadata: 
   Length     Class      Mode 
      141 character character 

}}
\keyword{datasets}
