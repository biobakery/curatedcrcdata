\name{ TCGA.RNASeqV2.READ_eset }
\alias{ TCGA.RNASeqV2.READ_eset }
\docType{data}
\title{ Comprehensive molecular characterization of human colon and rectal cancer. }
\description{ To characterize somatic alterations in colorectal carcinoma, we conducted a genome-scale analysis of 276 samples, analysing exome sequence, DNA copy number, promoter methylation and messenger RNA and microRNA expression. A subset of these samples (97) underwent low-depth-of-coverage whole-genome sequencing. In total, 16\% of colorectal carcinomas were found to be hypermutated: three-quarters of these had the expected high microsatellite instability, usually with hypermethylation and MLH1 silencing, and one-quarter had somatic mismatch-repair gene and polymerase ?? (POLE) mutations. Excluding the hypermutated cancers, colon and rectum cancers were found to have considerably similar patterns of genomic alteration. Twenty-four genes were significantly mutated, and in addition to the expected APC, TP53, SMAD4, PIK3CA and KRAS mutations, we found frequent mutations in ARID1A, SOX9 and FAM123B. Recurrent copy-number alterations include potentially drug-targetable amplifications of ERBB2 and newly discovered amplification of IGF2. Recurrent chromosomal translocations include the fusion of NAV2 and WNT pathway member TCF7L1. Integrative analyses suggest new markers for aggressive colorectal carcinoma and an important role for MYC-directed transcriptional activation and repression. }
\usage{data( TCGA.RNASeqV2.READ_eset )}
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
assayData: 20502 features, 6 samples
Platform type: NA
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   3 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
   3.00    3.00    3.00    3.00    3.44    2.72      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

unique_patient_ID: 
   Length     Class      Mode 
        6 character character 

alt_sample_name: 
   Length     Class      Mode 
        6 character character 

sample_type: 
tumor 
    6 

primarysite: 
re 
 6 

age_at_initial_pathologic_diagnosis: 
56 57 72 73 77 
 1  1  1  1  2 

days_to_death: 
 991 1257 1741 NA's 
   1    1    1    3 

vital_status: 
deceased   living     NA's 
       3        2        1 

msi: 
   n NA's 
   5    1 

location: 
Rectosigmoid Junction                rectum 
                    4                     2 

gender: 
f 
6 

kras: 
n y 
3 3 

mss: 
   y NA's 
   5    1 

uncurated_author_metadata: 
   Length     Class      Mode 
        6 character character 

}}
\keyword{datasets}
