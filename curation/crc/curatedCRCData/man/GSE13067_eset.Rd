\name{ GSE13067_eset }
\alias{ GSE13067_eset }
\docType{data}
\title{ DNA copy-number alterations underlie gene expression differences between microsatellite stable and unstable colorectal cancers. }
\description{ About 15\% of colorectal cancers harbor microsatellite instability (MSI). MSI-associated gene expression changes have been identified in colorectal cancers, but little overlap exists between signatures hindering an assessment of overall consistency. Little is known about the causes and downstream effects of differential gene expression.DNA microarray data on 89 MSI and 140 microsatellite-stable (MSS) colorectal cancers from this study and 58 MSI and 77 MSS cases from three published reports were randomly divided into test and training sets. MSI-associated gene expression changes were assessed for cross-study consistency using training samples and validated as MSI classifier using test samples. Differences in biological pathways were identified by functional category analysis. Causation of differential gene expression was investigated by comparison to DNA copy-number data.MSI-associated gene expression changes in colorectal cancers were found to be highly consistent across multiple studies of primary tumors and cancer cell lines from patients of different ethnicities (P < 0.001). Clustering based on consistent changes separated additional test cases by MSI status, and classification of individual samples predicted MSI status with a sensitivity of 96\% and specificity of 85\%. Genes associated with immune response were up-regulated in MSI cancers, whereas genes associated with cell-cell adhesion, ion binding, and regulation of metabolism were down-regulated. Differential gene expression was shown to reflect systematic differences in DNA copy-number aberrations between MSI and MSS tumors (P < 0.001).Our results show cross-study consistency of MSI-associated gene expression changes in colorectal cancers. DNA copy-number alterations partly cause the differences in gene expression between MSI and MSS cancers. }
\usage{data( GSE13067_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Jorissen RN, Lipton L, Gibbs P, Chapman M et al.??DNA copy-number alterations underlie gene expression differences between microsatellite stable and unstable colorectal cancers.??Clin Cancer Res??2008 Dec 15 
  Laboratory: Jorissen and Sieber 2008 
  Contact information:  
  Title: DNA copy-number alterations underlie gene expression differences between microsatellite stable and unstable colorectal cancers. 
  URL:  
  PMIDs: 19088021 

  Abstract: A 251 word abstract is available. Use 'abstract' method.
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
assayData: 19320 features, 74 samples
Platform type: hgu133plus2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       74 character character 

sample_type: 
tumor 
   74 

msi: 
MSI MSS 
 11  63 

batch: 
   Length     Class      Mode 
       74 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
       74 character character 

}}
\keyword{datasets}
