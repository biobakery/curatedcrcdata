\name{ GSE21510_eset }
\alias{ GSE21510_eset }
\docType{data}
\title{ Clinical significance of osteoprotegerin expression in human colorectal cancer. }
\description{ This study aimed to identify a novel biomarker or a target of treatment for colorectal cancer (CRC).The expression profiles of cancer cells in 104 patients with CRC were examined using laser microdissection and oligonucleotide microarray analysis. Overexpression in CRC cells, especially in patients with distant metastases, was a prerequisite to select candidate genes. The mRNA expression of candidate genes was investigated by quantitative reverse transcriptase PCR (RT-PCR) in 77 patients as a validation study. We analyzed the protein expression and localization of the candidate gene by immunohistochemical study and investigated the relationship between protein expression and clinicopathologic features in 274 CRC patients.Using microarray analysis, we identified 6 candidate genes related to distant metastases in CRC patients. Among these genes, osteoprotegerin (OPG) is known to be associated with aggressiveness in several cancers through inhibition of apoptosis via neutralization of the function of TNF-related apoptosis-inducing ligand. The mRNA expression of OPG in cancer tissues was significantly higher in patients with distant metastases than those without metastases. Overexpression of OPG protein was associated with significantly worse overall survival and relapse-free survival. Moreover, overexpression of the OPG protein was an independent risk factor for CRC recurrence.Overexpression of OPG may be a predictive biomarker of CRC recurrence and a target for treatment of this disease.??2011 AACR. }
\usage{data( GSE21510_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Tsukamoto S, Ishikawa T, Iida S, Ishiguro M et al.??Clinical significance of osteoprotegerin expression in human colorectal cancer.??Clin Cancer Res2011 Apr 15 
  Laboratory: Tsukamoto, Sugihara 2010 
  Contact information:  
  Title: Clinical significance of osteoprotegerin expression in human colorectal cancer. 
  URL:  
  PMIDs: 21270110 

  Abstract: A 211 word abstract is available. Use 'abstract' method.
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
assayData: 19320 features, 148 samples
Platform type: hgu133plus2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      148 character character 

sample_type: 
adjacentnormal          tumor 
            25            123 

summarystage: 
early  late  NA's 
   19    74    55 

M: 
  0   1 
121  27 

Dstage: 
   B    C    D NA's 
  54   47   27   20 

stageall: 
 0  1  2  3  4 
 1 19 54 47 27 

ethnicity: 
other 
  148 

batch: 
   Length     Class      Mode 
      148 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
      148 character character 

}}
\keyword{datasets}
