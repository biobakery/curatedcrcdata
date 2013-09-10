\name{ TCGA.RNASeqV2_eset }
\alias{ TCGA.RNASeqV2_eset }
\docType{data}
\title{ Comprehensive molecular characterization of human colon and rectal cancer. }
\description{ To characterize somatic alterations in colorectal carcinoma, we conducted a genome-scale analysis of 276 samples, analysing exome sequence, DNA copy number, promoter methylation and messenger RNA and microRNA expression. A subset of these samples (97) underwent low-depth-of-coverage whole-genome sequencing. In total, 16\% of colorectal carcinomas were found to be hypermutated: three-quarters of these had the expected high microsatellite instability, usually with hypermethylation and MLH1 silencing, and one-quarter had somatic mismatch-repair gene and polymerase ?? (POLE) mutations. Excluding the hypermutated cancers, colon and rectum cancers were found to have considerably similar patterns of genomic alteration. Twenty-four genes were significantly mutated, and in addition to the expected APC, TP53, SMAD4, PIK3CA and KRAS mutations, we found frequent mutations in ARID1A, SOX9 and FAM123B. Recurrent copy-number alterations include potentially drug-targetable amplifications of ERBB2 and newly discovered amplification of IGF2. Recurrent chromosomal translocations include the fusion of NAV2 and WNT pathway member TCF7L1. Integrative analyses suggest new markers for aggressive colorectal carcinoma and an important role for MYC-directed transcriptional activation and repression. }
\usage{data( TCGA.RNASeqV2_eset )}
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
      [RNASeqV2] Illumina HiSeq RNA sequencing
   platform_shorttitle:     
      
   platform_summary:     
      NA
   platform_manufacturer:     
      Illumina
   platform_distribution:     
      sequencing
   platform_accession:     
      NA
   platform_technology:     
      in situ oligonucleotide
   warnings:     
      No warnings yet

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: ? A1BG ... ZZZ3 (20502 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 20502 features, 195 samples
Platform type: NA
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   174 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 21.000  21.000  21.000  18.000   1.208   0.715   5.605 

--------------------------- 
Available sample meta-data: 
--------------------------- 

unique_patient_ID: 
   Length     Class      Mode 
      195 character character 

alt_sample_name: 
   Length     Class      Mode 
      195 character character 

sample_type: 
adjacentnormal          tumor 
            14            181 

primarysite: 
 co 
195 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  31.00   55.00   66.00   64.95   77.00   90.00 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   43.0   187.0   424.0   816.6  1422.0  2134.0     174 

vital_status: 
deceased   living     NA's 
      18      119       58 

msi: 
   n    y NA's 
  32    9  154 

location: 
   Length     Class      Mode 
      195 character character 

gender: 
  f   m 
 87 108 

kras: 
   n    y NA's 
   9    4  182 

mss: 
   n    y NA's 
   9   32  154 

uncurated_author_metadata: 
   Length     Class      Mode 
      195 character character 

fu: 
   n NA's 
  60  135 

bevacizumab: 
   n NA's 
  60  135 

irinotecan: 
   n NA's 
  60  135 

capecitabine: 
   n    y NA's 
  59    1  135 

cpt11: 
   n    y NA's 
  59    1  135 

dexamethasone: 
   n    y NA's 
  59    1  135 

erbitux: 
   n NA's 
  60  135 

gcsf: 
   n NA's 
  60  135 

fudr: 
   n    y NA's 
  57    3  135 

folfiri: 
   n    y NA's 
  59    1  135 

folfox: 
   n NA's 
  60  135 

leucovorin: 
   n NA's 
  60  135 

mitomycin: 
   n NA's 
  60  135 

platin: 
   n NA's 
  60  135 

panitumumab: 
   n NA's 
  60  135 

pegfilgrastim: 
   n    y NA's 
  57    3  135 

raltitrexed: 
   n    y NA's 
  59    1  135 

xeloda: 
   n NA's 
  60  135 

ancillary: 
   n    y NA's 
  58    2  135 

chemotherapy: 
   n    y NA's 
   2   58  135 

moltherapy: 
   n    y NA's 
  57    3  135 

}}
\keyword{datasets}
