\name{ GSE28702_eset }
\alias{ GSE28702_eset }
\docType{data}
\title{ Potential responders to FOLFOX therapy for colorectal cancer by Random Forests analysis. }
\description{ Molecular characterisation using gene-expression profiling will undoubtedly improve the prediction of treatment responses, and ultimately, the clinical outcome of cancer patients.To establish the procedures to identify responders to FOLFOX therapy, 83 colorectal cancer (CRC) patients including 42 responders and 41 non-responders were divided into training (54 patients) and test (29 patients) sets. Using Random Forests (RF) algorithm in the training set, predictor genes for FOLFOX therapy were identified, which were applied to test samples and sensitivity, specificity, and out-of-bag classification accuracy were calculated.In the training set, 22 of 27 responders (81.4\% sensitivity) and 23 of 27 non-responders (85.1\% specificity) were correctly classified. To improve the prediction model, we removed the outliers determined by RF, and the model could correctly classify 21 of 23 responders (91.3\%) and 22 of 23 non-responders (95.6\%) in the training set, and 80.0\% sensitivity and 92.8\% specificity, with an accuracy of 69.2\% in 29 independent test samples.Random Forests on gene-expression data for CRC patients was effectively able to stratify responders to FOLFOX therapy with high accuracy, and use of pharmacogenomics in anticancer therapy is the first step in planning personalised therapy. }
\usage{data( GSE28702_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Tsuji S, Midorikawa Y, Takahashi T, Yagi K et al.??Potential responders to FOLFOX therapy for colorectal cancer by Random Forests analysis.??Br J Cancer??2012 Jan 3 
  Laboratory: Tsuji 2011 
  Contact information:  
  Title: Potential responders to FOLFOX therapy for colorectal cancer by Random Forests analysis. 
  URL:  
  PMIDs: 22095227 

  Abstract: A 185 word abstract is available. Use 'abstract' method.
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
assayData: 19320 features, 59 samples
Platform type: hg-u133_plus_2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       59 character character 

sample_type: 
metastatic      tumor 
         3         56 

gender: 
 f  m 
20 39 

drug_name: 
mfolfox6 
      59 

drug_response: 
 n  y 
25 34 

uncurated_author_metadata: 
   Length     Class      Mode 
       59 character character 

}}
\keyword{datasets}
