\name{ GSE3294_eset }
\alias{ GSE3294_eset }
\docType{data}
\title{ Comparative study of gene expression by cDNA microarray in human colorectal cancer tissues and normal mucosa. }
\description{ The causative molecular pathways underlying the pathogenesis of colorectal cancer (CRC) need to be better characterized. The purpose of our study was to better understand the genetic mechanism of oncogenesis for human colorectal cancer and to identify new potential tumor markers of use in clinical practice. We used cDNA microarrays to compare gene expression profiles of colorectal biopsies from 25 CRC patients and 13 normal mucosa from adjacent non-cancerous tissues. Findings were validated by real-time PCR; in addition, western blotting and immunochemistry analysis were carried out as further confirmation of differential expression at a protein level. Comparing cancerous tissues with normal colonic mucosa we identified 584 known genes differentially expressed to a significant degree (p<0.001). Many of the transcripts that were more abundant in tumors than in non-neoplastic tissues appear to reflect important events for colon carcinogenesis. For example, a significant number of these genes serve as apoptotic inhibitors (e.g. BFAR, BIRC1, BIRC6). Furthermore, we observed the simultaneous up-regulation of HLA-E and the down-regulation of beta2-microglobulin; these genes strongly support a potential tumor escape strategy from immune surveillance in colon cancer tissues. Our study provides new gene candidates in the pathogenesis of human CRC disease. From our results we hypothesize that CRC cells escape immune surveillance through a specific gene expression alteration; moreover, over-expression of several survival genes seems to confer a more anti-apoptotic phenotype. These genes are involved in pathways not previously implicated in CRC pathogenesis and they may provide new targets for therapy. }
\usage{data( GSE3294_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Bianchini M, Levy E, Zucchini C, Pinski V et al.??Comparative study of gene expression by cDNA microarray in human colorectal cancer tissues and normal mucosa.??Int J Oncol??2006 Jul 
  Laboratory: Bianchini 2005 
  Contact information:  
  Title: Comparative study of gene expression by cDNA microarray in human colorectal cancer tissues and normal mucosa. 
  URL:  
  PMIDs: 16773188 

  Abstract: A 245 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      UHN SS-Human 19Kv7
   platform_shorttitle:     
      UHN SS-Human 19Kv7
   platform_summary:     
      uhn ss-human 19kv7
   platform_manufacturer:     
      UHN
   platform_distribution:     
      commercial
   platform_accession:     
      GPL2829
   platform_technology:     
      spotted DNA/cDNA
   warnings:     
      No warnings yet

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: AA001103 AA001104 ... Z45302 (15437 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 15437 features, 24 samples
Platform type: uhn ss-human 19kv7
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       24 character character 

summarystage: 
early  late 
    4    20 

T: 
 1  2  3 
 3  1 20 

N: 
 0  1  2 
10  9  5 

M: 
 0 
24 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  35.00   67.25   76.50   71.46   78.75   89.00 

gender: 
 f  m 
 8 16 

stageall: 
 2  3 
 6 18 

uncurated_author_metadata: 
   Length     Class      Mode 
       24 character character 

}}
\keyword{datasets}
