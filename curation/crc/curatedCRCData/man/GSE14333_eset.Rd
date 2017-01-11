\name{ GSE14333_eset }
\alias{ GSE14333_eset }
\docType{data}
\title{ Metastasis-Associated Gene Expression Changes Predict Poor Outcomes in Patients with Dukes Stage B and C Colorectal Cancer. }
\description{ PURPOSE: Colorectal cancer prognosis is currently predicted from pathologic staging, providing limited discrimination for Dukes stage B and C disease. Additional markers for outcome are required to help guide therapy selection for individual patients. EXPERIMENTAL DESIGN: A multisite single-platform microarray study was done on 553 colorectal cancers. Gene expression changes were identified between stage A and D tumors (three training sets) and assessed as a prognosis signature in stage B and C tumors (independent test and external validation sets). RESULTS: One hundred twenty-eight genes showed reproducible expression changes between three sets of stage A and D cancers. Using consistent genes, stage B and C cancers clustered into two groups resembling early-stage and metastatic tumors. A Prediction Analysis of Microarray algorithm was developed to classify individual intermediate-stage cancers into stage A-like/good prognosis or stage D-like/poor prognosis types. For stage B patients, the treatment adjusted hazard ratio for 6-year recurrence in individuals with stage D-like cancers was 10.3 (95\% confidence interval, 1.3-80.0; P = 0.011). For stage C patients, the adjusted hazard ratio was 2.9 (95\% confidence interval, 1.1-7.6; P = 0.016). Similar results were obtained for an external set of stage B and C patients. The prognosis signature was enriched for downregulated immune response genes and upregulated cell signaling and extracellular matrix genes. Accordingly, sparse tumor infiltration with mononuclear chronic inflammatory cells was associated with poor outcome in independent patients. CONCLUSIONS: Metastasis-associated gene expression changes can be used to refine traditional outcome prediction, providing a rational approach for tailoring treatments to subsets of patients. (Clin Cancer Res 2009;15(24):7642-51). }
\usage{data( GSE14333_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Jorissen RN, Gibbs P, Christie M, Prakash S et al.??Metastasis-Associated Gene Expression Changes Predict Poor Outcomes in Patients with Dukes Stage B and C Colorectal Cancer.??Clin Cancer Res??2009 Dec 15 
  Laboratory: Jorissen and Sieber 2008 
  Contact information:  
  Title: Metastasis-Associated Gene Expression Changes Predict Poor Outcomes in Patients with Dukes Stage B and C Colorectal Cancer. 
  URL:  
  PMIDs: 19996206 

  Abstract: A 257 word abstract is available. Use 'abstract' method.
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
assayData: 19320 features, 290 samples
Platform type: hgu133plus2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      290 character character 

sample_type: 
tumor 
  290 

primarysite: 
  co   re NA's 
 250   39    1 

summarystage: 
early  late  NA's 
   44   152    94 

N: 
   0 NA's 
 138  152 

M: 
  0   1 
229  61 

Dstage: 
 A  B  C  D 
44 94 91 61 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  26.00   58.00   67.00   65.96   75.00   92.00 

summarylocation: 
   l    r NA's 
 161  125    4 

gender: 
  f   m 
126 164 

stageall: 
   1    3    4 NA's 
  44   91   61   94 

dfs_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    50                    176                     64 

days_to_recurrence_or_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   27.6   668.5  1154.0  1306.0  1785.0  4276.0      64 

batch: 
   Length     Class      Mode 
      290 character character 

drug_treatment: 
   n    y NA's 
 172  117    1 

preop_drug_treatment: 
  n 
290 

chemotherapy: 
   n    y NA's 
 172  117    1 

uncurated_author_metadata: 
   Length     Class      Mode 
      290 character character 

}}
\keyword{datasets}
