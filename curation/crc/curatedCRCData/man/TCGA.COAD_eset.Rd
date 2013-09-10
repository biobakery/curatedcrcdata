\name{ TCGA.COAD_eset }
\alias{ TCGA.COAD_eset }
\docType{data}
\title{ Comprehensive molecular characterization of human colon and rectal cancer. }
\description{ To characterize somatic alterations in colorectal carcinoma, we conducted a genome-scale analysis of 276 samples, analysing exome sequence, DNA copy number, promoter methylation and messenger RNA and microRNA expression. A subset of these samples (97) underwent low-depth-of-coverage whole-genome sequencing. In total, 16\% of colorectal carcinomas were found to be hypermutated: three-quarters of these had the expected high microsatellite instability, usually with hypermethylation and MLH1 silencing, and one-quarter had somatic mismatch-repair gene and polymerase ?? (POLE) mutations. Excluding the hypermutated cancers, colon and rectum cancers were found to have considerably similar patterns of genomic alteration. Twenty-four genes were significantly mutated, and in addition to the expected APC, TP53, SMAD4, PIK3CA and KRAS mutations, we found frequent mutations in ARID1A, SOX9 and FAM123B. Recurrent copy-number alterations include potentially drug-targetable amplifications of ERBB2 and newly discovered amplification of IGF2. Recurrent chromosomal translocations include the fusion of NAV2 and WNT pathway member TCF7L1. Integrative analyses suggest new markers for aggressive colorectal carcinoma and an important role for MYC-directed transcriptional activation and repression. }
\usage{data( TCGA.COAD_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Comprehensive molecular characterization of human colon and rectal cancer. Nature 2012, 487:330-337 
  Laboratory: The Cancer Genome Atlas Network 2012 
  Contact information:  
  Title: Comprehensive molecular characterization of human colon and rectal cancer. 
  URL:  
  PMIDs: 22810696 

  Abstract: A 168 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Agilent 244K Custom Gene Expression G4502A-07-3
   platform_shorttitle:     
      Agilent G4502A-07-3
   platform_summary:     
      agilent-014850 whole human genome microarray 4x44k g4112f
   platform_manufacturer:     
      Agilent
   platform_distribution:     
      commercial
   platform_accession:     
      NA
   platform_technology:     
      in situ oligonucleotide
   warnings:     
      No warnings yet

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 15E1.2 2'-PDE ... ZZZ3 (17814 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 17814 features, 130 samples
Platform type: agilent-014850 whole human genome microarray 4x44k g4112f
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   124 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
   6.00    6.00    6.00    5.00    1.37    1.16      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

unique_patient_ID: 
   Length     Class      Mode 
      130 character character 

alt_sample_name: 
   Length     Class      Mode 
      130 character character 

sample_type: 
adjacentnormal          tumor 
             7            123 

primarysite: 
  co NA's 
 129    1 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  36.00   65.00   72.00   70.89   80.00   90.00 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   29.0   442.0   500.0   494.5   681.0   774.0     124 

vital_status: 
deceased   living     NA's 
       5       15      110 

msi: 
   n    y NA's 
   4    1  125 

location: 
   Length     Class      Mode 
      130 character character 

gender: 
 f  m 
67 63 

kras: 
   n    y NA's 
   1    2  127 

mss: 
   n    y NA's 
   1    4  125 

uncurated_author_metadata: 
   Length     Class      Mode 
      130 character character 

fu: 
   n NA's 
  31   99 

bevacizumab: 
   n NA's 
  31   99 

irinotecan: 
   n NA's 
  31   99 

capecitabine: 
   n    y NA's 
  25    6   99 

cpt11: 
   n NA's 
  31   99 

dexamethasone: 
   n NA's 
  31   99 

erbitux: 
   n NA's 
  31   99 

gcsf: 
   n NA's 
  31   99 

fudr: 
   n NA's 
  31   99 

folfiri: 
   n NA's 
  31   99 

folfox: 
   n NA's 
  31   99 

leucovorin: 
   n NA's 
  31   99 

mitomycin: 
   n    y NA's 
  30    1   99 

platin: 
   n NA's 
  31   99 

panitumumab: 
   n NA's 
  31   99 

pegfilgrastim: 
   n NA's 
  31   99 

raltitrexed: 
   n NA's 
  31   99 

xeloda: 
   n NA's 
  31   99 

ancillary: 
   n    y NA's 
  30    1   99 

chemotherapy: 
   y NA's 
  31   99 

moltherapy: 
   n    y NA's 
  21   10   99 

}}
\keyword{datasets}
