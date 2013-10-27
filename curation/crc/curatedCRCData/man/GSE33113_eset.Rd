\name{ GSE33113_eset }
\alias{ GSE33113_eset }
\docType{data}
\title{ Mutations in the Ras-Raf Axis underlie the prognostic value of CD133 in colorectal cancer. }
\description{ High expression of cancer stem cell (CSC) marker CD133 has been used as a predictor for prognosis in colorectal cancer (CRC), suggesting that enumeration of CSCs, using CD133, is predictive for disease progression. However, we showed recently that both CD133 mRNA and protein are not downregulated during differentiation of colon CSCs, pointing to an alternative reason for the prognostic value of CD133. We therefore set out to delineate the relation between CD133 expression and prognosis.A CRC patient series was studied for expression of CD133 and other CSC markers by microarray and quantitative PCR analysis. In addition, several common mutations were analyzed to determine the relation with CD133 expression.CD133 mRNA expression predicted relapse-free survival in our patient series, whereas several other CSC markers could not. Moreover, no correlation was found between expression of other CSC markers and CD133. Interestingly, high CD133 expression was related to mutations in K-Ras and B-Raf, and inhibition of mutant K-Ras or downstream mitogen-activated protein kinase kinase (MEK) signaling decreases CD133 expression. In addition, an activated K-Ras gene expression signature could predict CD133 expression in our patient set as well as data sets of other tumor types.CD133 expression is upregulated in CRC tumors that have a hyperactivated Ras-Raf-MEK-ERK pathway and is therefore related to mutations in K-Ras or B-Raf. As mutations in either gene have been related to poor prognosis, we conclude that CD133 expression is not indicative for CSC numbers but rather related to the mutation or activity status of the Ras-Raf pathway. }
\usage{data( GSE33113_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: de Sousa E Melo F, Colak S, Buikhuisen J, Koster J et al.??Methylation of cancer-stem-cell-associated Wnt target genes predicts poor prognosis in colorectal cancer patients.??Cell Stem Cell??2011 Nov 4 
  Laboratory: Medema JP,??Tanis PJ 2011 
  Contact information:  
  Title: Mutations in the Ras-Raf Axis underlie the prognostic value of CD133 in colorectal cancer. 
  URL:  
  PMIDs: 22496204 

  Abstract: A 247 word abstract is available. Use 'abstract' method.
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
assayData: 19320 features, 96 samples
Platform type: hg-u133_plus_2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       96 character character 

sample_type: 
adjacentnormal          tumor 
             6             90 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  35.00   60.00   72.50   69.73   79.00   95.00       2 

gender: 
   f    m NA's 
  50   44    2 

stageall: 
 2 
96 

days_to_recurrence_or_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   50.0   508.8  1180.0  1234.0  1633.0  3599.0       2 

batch: 
   Length     Class      Mode 
       96 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
       96 character character 

}}
\keyword{datasets}
