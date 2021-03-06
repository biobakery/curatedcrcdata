\name{ GSE26682.GPL96_eset }
\alias{ GSE26682.GPL96_eset }
\docType{data}
\title{ MRE11 deficiency increases sensitivity to poly(ADP-ribose) polymerase inhibition in microsatellite unstable colorectal cancers. }
\description{ Microsatellite instability (MSI) is displayed by approximately 15\% of colorectal cancers (CRC). Defective DNA mismatch repair generates mutations at repetitive DNA sequences such as those located in the double strand break (DSB) repair gene MRE11. We assessed the mutational status of MRE11 in a panel of 17 CRC cell lines and 46 primary tumors and found a strong correlation with MSI status in both cell lines and tumors. Therefore, we hypothesized that deficiency in MRE11 may sensitize CRC cells to poly(ADP-ribose) polymerase (PARP-1) inhibition based on the concept of synthetic lethality. We further assessed the activity of the PARP-1 inhibitor, ABT-888, in CRC cell lines and observed preferential cytotoxicity in those MSI cell lines harboring mutations in MRE11 compared with both wild-type cell lines and microsatellite stable (MSS) cell lines. A significant correlation between MRE11 expression levels and cytotoxicity to ABT-888 at 10 ??M was observed (R?? = 0.915, P < 0.001). Using two experimental approaches, including short hairpin RNA knocking down MRE11 in the wild-type and MSS cell line SW-480 and a second cell line model transfected with mutant MRE11, we experimentally tried to confirm the role of MRE11 in conferring sensitivity to PARP-1 inhibition. Both models led to changes in proliferation in response to ABT-888 at different concentrations, and a drug-response effect was not observed, suggesting a possible contribution of additional genes. We conclude that MSI colorectal tumors deficient in DSB repair secondary to mutation in MRE11 show a higher sensitivity to PARP-1 inhibition. Further clinical investigation of PARP-1 inhibitors is warranted in MSI CRCs. }
\usage{data( GSE26682.GPL96_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Vilar E,??Gruber SB,??Rennert G,??Bartnik CM,??Stenzel SL,??Iniesta MD,??Raskin L,Mukherjee B,??Ahn J,??Moreno V,??Morgan MA. MRE11 deficiency increases sensitivity to poly(ADP-ribose) polymerase inhibition in microsatellite unstable colorectal cancers. Cancer Res 2011 Apr 1 
  Laboratory: Vilar E,??Morgan MA 2011 
  Contact information:  
  Title: MRE11 deficiency increases sensitivity to poly(ADP-ribose) polymerase inhibition in microsatellite unstable colorectal cancers. 
  URL:  
  PMIDs: 21300766 

  Abstract: A 257 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG-U133A] Affymetrix Human Genome U133A Array
   platform_shorttitle:     
      Affymetrix HG-U133A
   platform_summary:     
      hgu133a
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL96
   platform_technology:     
      in situ oligonucleotide
   warnings:     
      No warnings yet

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (12986 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 12986 features, 155 samples
Platform type: hgu133a
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      155 character character 

sample_type: 
tumor 
  155 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  21.00   66.50   75.00   72.61   80.50   94.00 

msi: 
 MSI  MSS NA's 
  17  123   15 

gender: 
 f  m 
69 86 

batch: 
   Length     Class      Mode 
      155 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
      155 character character 

}}
\keyword{datasets}
