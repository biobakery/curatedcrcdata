\name{ GSE14095_eset }
\alias{ GSE14095_eset }
\docType{data}
\title{ Gene expression signature and response to the use of leucovorin, fluorouracil and oxaliplatin in colorectal cancer patients. }
\description{ FOLFOX (a combination of leucovorin, fluorouracil and oxaliplatin) has achieved substantial success in the treatment of colorectal cancer (CRC) patients. However, about half of all patients show resistance to this regimen and some develop adverse symptoms such as neurotoxicity. In order to select patients who would benefit most from this therapy, we aimed to build a predictor for the response to FOLFOX using microarray gene expression profiles of primary CRC samples.Forty patients who underwent surgery for primary lesions were examined. All patients had metastatic or recurrent CRC and received modified FOLFOX6. Responders and nonresponders were determined according to the best observed response at the end of the first-line treatment. Gene-expression profiles of primary CRC were determined using Human Genome GeneChip arrays U133. We identified discriminating genes whose expression differed significantly between responders and nonresponders and then carried out supervised class prediction using the k-nearest-neighbour method.We identified 27 probes that were differentially expressed between responders and nonresponders at significant levels. Based on the expression of these genes, we constructed a FOLFOX response predictor with an overall accuracy of 92.5\%. The sensitivity, specificity, positive and negative predictive values were 78.6\%, 100\%, 100\% and 89.7\%, respectively.The present model suggests the possibility of selecting patients who would benefit from FOLFOX therapy both in the metastatic and the adjuvant setting. To our knowledge, this is the first study to establish a prediction model for the response to FOLFOX chemotherapy based on gene expression by microarray analysis. }
\usage{data( GSE14095_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Watanabe T, Kobunai T, Yamamoto Y, Matsuda K et al.??Gene expression signature and response to the use of leucovorin, fluorouracil and oxaliplatin in colorectal cancer patients.??Clin Transl Oncol??2011 Jun 
  Laboratory: Watanabe, Hashimoto 2008 
  Contact information:  
  Title: Gene expression signature and response to the use of leucovorin, fluorouracil and oxaliplatin in colorectal cancer patients. 
  URL:  
  PMIDs: 21680303 

  Abstract: A 241 word abstract is available. Use 'abstract' method.
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
assayData: 19320 features, 189 samples
Platform type: hg-u133_plus_2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      189 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
      189 character character 

}}
\keyword{datasets}
