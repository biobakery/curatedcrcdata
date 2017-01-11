\name{ TCGA.READ_eset }
\alias{ TCGA.READ_eset }
\docType{data}
\title{ Comprehensive molecular characterization of human colon and rectal cancer. }
\description{ To characterize somatic alterations in colorectal carcinoma, we conducted a genome-scale analysis of 276 samples, analysing exome sequence, DNA copy number, promoter methylation and messenger RNA and microRNA expression. A subset of these samples (97) underwent low-depth-of-coverage whole-genome sequencing. In total, 16\% of colorectal carcinomas were found to be hypermutated: three-quarters of these had the expected high microsatellite instability, usually with hypermethylation and MLH1 silencing, and one-quarter had somatic mismatch-repair gene and polymerase ?? (POLE) mutations. Excluding the hypermutated cancers, colon and rectum cancers were found to have considerably similar patterns of genomic alteration. Twenty-four genes were significantly mutated, and in addition to the expected APC, TP53, SMAD4, PIK3CA and KRAS mutations, we found frequent mutations in ARID1A, SOX9 and FAM123B. Recurrent copy-number alterations include potentially drug-targetable amplifications of ERBB2 and newly discovered amplification of IGF2. Recurrent chromosomal translocations include the fusion of NAV2 and WNT pathway member TCF7L1. Integrative analyses suggest new markers for aggressive colorectal carcinoma and an important role for MYC-directed transcriptional activation and repression. }
\usage{data( TCGA.READ_eset )}
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
assayData: 17814 features, 51 samples
Platform type: agilent-014850 whole human genome microarray 4x44k g4112f
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   50 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
  1.000   1.000   1.000   1.000   0.866      NA      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

unique_patient_ID: 
   Length     Class      Mode 
       51 character character 

alt_sample_name: 
   Length     Class      Mode 
       51 character character 

sample_type: 
tumor 
   51 

primarysite: 
  re NA's 
  50    1 

summarystage: 
early  late 
   30    21 

T: 
 1  2  3  4 
 3  9 37  2 

N: 
 0  1  2 
30 11 10 

M: 
 0  1 
43  8 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  41.00   62.00   67.00   67.67   73.50   89.00 

days_to_tumor_recurrence: 
 301  327 NA's 
   1    1   49 

recurrence_status: 
norecurrence   recurrence         NA's 
           3            2           46 

days_to_death: 
 316 NA's 
   1   50 

vital_status: 
deceased   living     NA's 
       1        4       46 

msi: 
 MSS NA's 
   1   50 

location: 
rectosigmoid       rectum         NA's 
           3           47            1 

summarylocation: 
   l NA's 
  50    1 

gender: 
 f  m 
24 27 

kras: 
mutant     wt   NA's 
     2      2     47 

lymphnodesremoved: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  12.00   17.00   22.00   24.73   32.00   65.00 

lymphnodesinvaded: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.000   0.000   0.000   2.608   3.000  31.000 

stageall: 
   1    2    3    4 NA's 
  12   18   13    7    1 

ethnicity: 
    black caucasian      NA's 
        1         4        46 

drug_treatment: 
   y NA's 
   2   49 

preop_drug_treatment: 
 n 
51 

fu: 
   n NA's 
  16   35 

bevacizumab: 
   n NA's 
  16   35 

irinotecan: 
   n NA's 
  16   35 

cetuximab: 
   n NA's 
  16   35 

capecitabine: 
   n NA's 
  16   35 

dexamethasone: 
   n NA's 
  16   35 

gcsf: 
   n NA's 
  16   35 

fudr: 
   n NA's 
  16   35 

folfiri: 
   n NA's 
  16   35 

folfox: 
   n NA's 
  16   35 

leucovorin: 
   n NA's 
  16   35 

mitomycin: 
   n NA's 
  16   35 

platin: 
   n NA's 
  16   35 

panitumumab: 
   n NA's 
  16   35 

pegfilgrastim: 
   n NA's 
  16   35 

raltitrexed: 
   n NA's 
  16   35 

ancillary: 
   n NA's 
  17   34 

chemotherapy: 
   y NA's 
  17   34 

moltherapy: 
   n    y NA's 
  15    2   34 

uncurated_author_metadata: 
   Length     Class      Mode 
       51 character character 

}}
\keyword{datasets}
