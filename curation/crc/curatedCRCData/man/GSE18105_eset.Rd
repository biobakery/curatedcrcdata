\name{ GSE18105_eset }
\alias{ GSE18105_eset }
\docType{data}
\title{ MUC12 mRNA expression is an independent marker of prognosis in stage II and stage III colorectal cancer. }
\description{ Distant metastasis is the major cause of death in colorectal cancer (CRC) patients. To identify genes influencing the prognosis of patients with CRC, we compared gene expression in primary tumors with and without distant metastasis using an oligonucleotide microarray. We also examined the expression of the candidate gene in 100 CRC patients by quantitative real-time reverse transcription PCR and studied the relationship between its expression and the prognosis of patients with CRC. As a result, we identified MUC12 as a candidate gene involved in metastasis processes by microarray analysis. Quantitative real-time reverse transcription PCR showed that MUC12 expression was significantly lower in cancer tissues than in adjacent normal tissues (p < 0.001). In Stages II and III CRC, patients with low expression showed worse disease-free survival (p = 0.020). Multivariate analysis disclosed that MUC12 expression status was an independent prognostic factor in Stages II and III CRC (relative risk, 8.236; 95\% confidence interval, 1.702-39.849 p = 0.009). Our study revealed the prognostic value of MUC12 expression in CRC patients. Moreover, our result suggests MUC12 expression is a possible candidate gene for assessing postoperative adjuvant therapy for CRC patients. }
\usage{data( GSE18105_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Matsuyama T, Ishikawa T, Mogushi K, Yoshida T et al.??MUC12 mRNA expression is an independent marker of prognosis in stage II and stage III colorectal cancer.??Int J Cancer??2010 Nov 15 
  Laboratory: Matsuyama, Sugihara 2009 
  Contact information:  
  Title: MUC12 mRNA expression is an independent marker of prognosis in stage II and stage III colorectal cancer. 
  URL:  
  PMIDs: 20162577 

  Abstract: A 188 word abstract is available. Use 'abstract' method.
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
assayData: 19320 features, 111 samples
Platform type: hg-u133_plus_2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      111 character character 

sample_type: 
adjacentnormal          tumor 
            17             94 

M: 
 0  1 
67 44 

uncurated_author_metadata: 
   Length     Class      Mode 
      111 character character 

}}
\keyword{datasets}
