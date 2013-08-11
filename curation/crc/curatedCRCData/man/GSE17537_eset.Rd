\name{ GSE17537_eset }
\alias{ GSE17537_eset }
\docType{data}
\title{ Experimentally derived metastasis gene expression profile predicts recurrence and death in patients with colon cancer. }
\description{ Staging inadequately predicts metastatic risk in patients with colon cancer. We used a gene expression profile derived from invasive, murine colon cancer cells that were highly metastatic in an immunocompetent mouse model to identify patients with colon cancer at risk of recurrence.This phase 1, exploratory biomarker study used 55 patients with colorectal cancer from Vanderbilt Medical Center (VMC) as the training dataset and 177 patients from the Moffitt Cancer Center as the independent dataset. The metastasis-associated gene expression profile developed from the mouse model was refined with comparative functional genomics in the VMC gene expression profiles to identify a 34-gene classifier associated with high risk of metastasis and death from colon cancer. A metastasis score derived from the biologically based classifier was tested in the Moffitt dataset.A high score was significantly associated with increased risk of metastasis and death from colon cancer across all pathologic stages and specifically in stage II and stage III patients. The metastasis score was shown to independently predict risk of cancer recurrence and death in univariate and multivariate models. For example, among stage III patients, a high score translated to increased relative risk of cancer recurrence (hazard ratio, 4.7; 95\% confidence interval, 1.566-14.05). Furthermore, the metastasis score identified patients with stage III disease whose 5-year recurrence-free survival was >88\% and for whom adjuvant chemotherapy did not increase survival time.A gene expression profile identified from an experimental model of colon cancer metastasis predicted cancer recurrence and death, independently of conventional measures, in patients with colon cancer.Copyright 2010 AGA Institute. Published by Elsevier Inc. All rights reserved. }
\usage{data( GSE17537_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information:  
  Title: Experimentally derived metastasis gene expression profile predicts recurrence and death in patients with colon cancer. 
  URL:  
  PMIDs: 19914252 

  Abstract: A 260 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
   platform_shorttitle:     
      
   platform_summary:     
      hg-u133_plus_2
   platform_manufacturer:     
      
   platform_distribution:     
      
   platform_accession:     
      
   platform_technology:     
      1

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1BG-AS1 ... ZZZ3 (19320 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19320 features, 55 samples
Platform type: hg-u133_plus_2
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
     55      55      55      20      NA      NA      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       55 character character 

summarygrade: 
high  low NA's 
   3   33   19 

G: 
   1    2    3 NA's 
   1   32    3   19 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  23.00   54.00   62.00   62.31   72.00   94.00 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.00   13.81 1243.00  988.10 1694.00 2303.00 

recurrence_status: 
norecurrence   recurrence 
          36           19 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  12.82  950.30 1506.00 1357.00 1801.00 3345.00 

vital_status: 
deceased   living 
      20       35 

gender: 
 f  m 
29 26 

stageall: 
 1  2  3  4 
 4 15 19 17 

ethnicity: 
    black caucasian  hispanic 
        4        50         1 

uncurated_author_metadata: 
   Length     Class      Mode 
       55 character character 

}}
\keyword{datasets}
