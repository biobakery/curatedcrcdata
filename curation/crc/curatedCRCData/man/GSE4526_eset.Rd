\name{ GSE4526_eset }
\alias{ GSE4526_eset }
\docType{data}
\title{ Gene expression signature for recurrence in stage III colorectal cancers. }
\description{ Colorectal cancer patients with lymph node metastases (stage III) show poorer prognosis than those without. Predicting development of recurrence may guide the need for intensive follow-up and/or adjuvant chemotherapy in such patients. The authors' objective was to identify a set of discriminating genes that could predict recurrence in stage III colorectal cancer.Thirty-six stage III colorectal cancer patients were studied. Tumor samples were obtained from surgically resected specimens. Thirteen patients developed recurrence, whereas 23 patients did not. Gene expression profiles were determined using human HG-U133 Plus 2.0 Gene Chip (Affymetrix, Santa Clara, Calif).The authors identified 45 discriminating genes between patients with and without recurrence. By using this gene set, they established a new model to predict recurrence with an accuracy of 90.9\%. The discriminating genes included calcineurin-binding protein 1 (CABIN1), whose expression differed remarkably between patients with and without recurrence (P=.0073). The authors further examined the DNA copy number of CABIN1 and were able to show a significant relation with recurrence (P<.012). Patients having CABIN1 gene loss demonstrated a higher risk of recurrence (odds ratio, 18.8). DNA copy number of CABIN1 alone could predict recurrence with an accuracy of 80.0\%.The results of the current study demonstrated that gene expression profiling is useful in predicting recurrence in stage III colorectal cancer. The authors identified CABIN1 among discriminating genes that may play a key role in the development of recurrence. These results may help to establish an individualized therapy for stage III colorectal cancer.Copyright (c) 2009 American Cancer Society. }
\usage{data( GSE4526_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Watanabe T, Kobunai T, Sakamoto E, Yamamoto Y et al.??Gene expression signature for recurrence in stage III colorectal cancers.??Cancer??2009 Jan 15 
  Laboratory: Watanabe T,??Kobunai T,??Toda E,??Oka T 2006 
  Contact information:  
  Title: Gene expression signature for recurrence in stage III colorectal cancers. 
  URL:  
  PMIDs: 19016304 

  Abstract: A 246 word abstract is available. Use 'abstract' method.
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

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1BG-AS1 ... ZZZ3 (19320 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19320 features, 36 samples
Platform type: hg-u133_plus_2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       36 character character 

recurrence_status: 
norecurrence   recurrence 
          23           13 

stageall: 
 3 
36 

uncurated_author_metadata: 
   Length     Class      Mode 
       36 character character 

}}
\keyword{datasets}
