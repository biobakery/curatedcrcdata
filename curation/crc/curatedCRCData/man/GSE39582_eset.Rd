\name{ GSE39582_eset }
\alias{ GSE39582_eset }
\docType{data}
\title{ Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. }
\description{ Colon cancer (CC) pathological staging fails to accurately predict recurrence, and to date, no gene expression signature has proven reliable for prognosis stratification in clinical practice, perhaps because CC is a heterogeneous disease. The aim of this study was to establish a comprehensive molecular classification of CC based on mRNA expression profile analyses.Fresh-frozen primary tumor samples from a large multicenter cohort of 750 patients with stage I to IV CC who underwent surgery between 1987 and 2007 in seven centers were characterized for common DNA alterations, including BRAF, KRAS, and TP53 mutations, CpG island methylator phenotype, mismatch repair status, and chromosomal instability status, and were screened with whole genome and transcriptome arrays. 566 samples fulfilled RNA quality requirements. Unsupervised consensus hierarchical clustering applied to gene expression data from a discovery subset of 443 CC samples identified six molecular subtypes. These subtypes were associated with distinct clinicopathological characteristics, molecular alterations, specific enrichments of supervised gene expression signatures (stem cell phenotype-like, normal-like, serrated CC phenotype-like), and deregulated signaling pathways. Based on their main biological characteristics, we distinguished a deficient mismatch repair subtype, a KRAS mutant subtype, a cancer stem cell subtype, and three chromosomal instability subtypes, including one associated with down-regulated immune pathways, one with up-regulation of the Wnt pathway, and one displaying a normal-like gene expression profile. The classification was validated in the remaining 123 samples plus an independent set of 1,058 CC samples, including eight public datasets. Furthermore, prognosis was analyzed in the subset of stage II-III CC samples. The subtypes C4 and C6, but not the subtypes C1, C2, C3, and C5, were independently associated with shorter relapse-free survival, even after adjusting for age, sex, stage, and the emerging prognostic classifier Oncotype DX Colon Cancer Assay recurrence score (hazard ratio 1.5, 95\% CI 1.1-2.1, p???=???0.0097). However, a limitation of this study is that information on tumor grade and number of nodes examined was not available.We describe the first, to our knowledge, robust transcriptome-based classification of CC that improves the current disease stratification based on clinicopathological variables and common DNA markers. The biological relevance of these subtypes is illustrated by significant differences in prognosis. This analysis provides possibilities for improving prognostic models and therapeutic strategies. In conclusion, we report a new classification of CC into six molecular subtypes that arise through distinct biological pathways. }
\usage{data( GSE39582_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Marisa L, de Reyni??A, Duval A, Selves J et al.??Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value.??PLoS Med??2013 
  Laboratory: Marisa, Boige 2012 
  Contact information:  
  Title: Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. 
  URL:  
  PMIDs: 23700391 

  Abstract: A 384 word abstract is available. Use 'abstract' method.
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
assayData: 19320 features, 566 samples
Platform type: hg-u133_plus_2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      566 character character 

T: 
  0   1   2   3   4 
  4  33 264 205  60 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   22.0    59.0    68.0    66.9    76.0    97.0       1 

location: 
  distal proximal 
     342      224 

gender: 
  f   m 
256 310 

stageall: 
  2 
566 

batch: 
   Length     Class      Mode 
      566 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
      566 character character 

}}
\keyword{datasets}
