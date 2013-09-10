\name{ GSE3964_eset }
\alias{ GSE3964_eset }
\docType{data}
\title{ Deciphering cellular states of innate tumor drug responses. }
\description{ The molecular mechanisms underlying innate tumor drug resistance, a major obstacle to successful cancer therapy, remain poorly understood. In colorectal cancer (CRC), molecular studies have focused on drug-selected tumor cell lines or individual candidate genes using samples derived from patients already treated with drugs, so that very little data are available prior to drug treatment.Transcriptional profiles of clinical samples collected from CRC patients prior to their exposure to a combined chemotherapy of folinic acid, 5-fluorouracil and irinotecan were established using microarrays. Vigilant experimental design, power simulations and robust statistics were used to restrain the rates of false negative and false positive hybridizations, allowing successful discrimination between drug resistance and sensitivity states with restricted sampling. A list of 679 genes was established that intrinsically differentiates, for the first time prior to drug exposure, subsequently diagnosed chemo-sensitive and resistant patients. Independent biological validation performed through quantitative PCR confirmed the expression pattern on two additional patients. Careful annotation of interconnected functional networks provided a unique representation of the cellular states underlying drug responses.Molecular interaction networks are described that provide a solid foundation on which to anchor working hypotheses about mechanisms underlying in vivo innate tumor drug responses. These broad-spectrum cellular signatures represent a starting point from which by-pass chemotherapy schemes, targeting simultaneously several of the molecular mechanisms involved, may be developed for critical therapeutic intervention in CRC patients. The demonstrated power of this research strategy makes it generally applicable to other physiological and pathological situations. }
\usage{data( GSE3964_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Graudens E, Boulanger V, Mollard C, Mariage-Samson R et al.??Deciphering cellular states of innate tumor drug responses.??Genome Biol??2006 
  Laboratory: Graudens, Imbeaud 2006 
  Contact information:  
  Title: Deciphering cellular states of innate tumor drug responses. 
  URL:  
  PMIDs: 16542501 

  Abstract: A 242 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      11K_VJF-ARRAY
   platform_shorttitle:     
      NA
   platform_summary:     
      11k_vjf-array
   platform_manufacturer:     
      3
   platform_distribution:     
      non-commercial
   platform_accession:     
      GPL3282
   platform_technology:     
      spotted DNA/cDNA
   warnings:     
      No warnings yet

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 384D8-2 38600 ... ZYX (5845 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 5845 features, 15 samples
Platform type: 11k_vjf-array
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       15 character character 

sample_type: 
adjacentnormal          tumor 
             6              9 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  49.00   60.50   65.00   64.27   68.00   71.00 

gender: 
f m 
7 8 

stageall: 
 4 
15 

uncurated_author_metadata: 
   Length     Class      Mode 
       15 character character 

}}
\keyword{datasets}
