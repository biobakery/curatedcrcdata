\name{ GSE12945_eset }
\alias{ GSE12945_eset }
\docType{data}
\title{ An expression module of WIPF1-coexpressed genes identifies patients with favorable prognosis in three tumor types. }
\description{ Wiskott-Aldrich syndrome (WAS) predisposes patients to leukemia and lymphoma. WAS is caused by mutations in the protein WASP which impair its interaction with the WIPF1 protein. Here, we aim to identify a module of WIPF1-coexpressed genes and to assess its use as a prognostic signature for colorectal cancer, glioma, and breast cancer patients. Two public colorectal cancer microarray data sets were used for discovery and validation of the WIPF1 co-expression module. Based on expression of the WIPF1 signature, we classified more than 400 additional tumors with microarray data from our own experiments or from publicly available data sets according to their WIPF1 signature expression. This allowed us to separate patient populations for colorectal cancers, breast cancers, and gliomas for which clinical characteristics like survival times and times to relapse were analyzed. Groups of colorectal cancer, breast cancer, and glioma patients with low expression of the WIPF1 co-expression module generally had a favorable prognosis. In addition, the majority of WIPF1 signature genes are individually correlated with disease outcome in different studies. Literature gene network analysis revealed that among WIPF1 co-expressed genes known direct transcriptional targets of c-myc, ESR1 and p53 are enriched. The mean expression profile of WIPF1 signature genes is correlated with the profile of a proliferation signature. The WIPF1 signature is the first microarray-based prognostic expression signature primarily developed for colorectal cancer that is instrumental in other tumor types: low expression of the WIPF1 module is associated with better prognosis. }
\usage{data( GSE12945_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Staub E, Groene J, Heinze M, Mennerich D et al.??An expression module of WIPF1-coexpressed genes identifies patients with favorable prognosis in three tumor types.??J Mol Med (Berl)??2009 Jun 
  Laboratory: Staub, Rosenthal 2009 
  Contact information:  
  Title: An expression module of WIPF1-coexpressed genes identifies patients with favorable prognosis in three tumor types. 
  URL:  
  PMIDs: 19399471 

  Abstract: A 241 word abstract is available. Use 'abstract' method.
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
assayData: 12986 features, 62 samples
Platform type: hgu133a
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
     62      62      62      12      NA      NA      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       62 character character 

sample_type: 
tumor 
   62 

primarysite: 
co re 
29 33 

summarygrade: 
high  low 
  31   31 

G: 
 2  3 
31 31 

summarystage: 
early  late 
   36    26 

T: 
 2  3  4 
16 42  4 

N: 
 0  1  2 
36 14 12 

M: 
 0  1  X 
56  5  1 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  38.00   59.00   65.00   64.45   73.75   87.00 

recurrence_status: 
norecurrence   recurrence         NA's 
          47            8            7 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    210    1005    1395    1267    1620    1920 

vital_status: 
deceased   living 
      12       50 

location: 
   Length     Class      Mode 
       62 character character 

summarylocation: 
   l    r NA's 
  48   13    1 

gender: 
 f  m 
28 34 

lymphnodesremoved: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  11.00   13.25   16.50   19.00   22.75   42.00 

lymphnodesinvaded: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.000   0.000   0.000   1.903   1.750  20.000 

stageall: 
 1  2  3  4 
13 23 21  5 

dfs_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    15                     46                      1 

days_to_recurrence_or_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     90    1110    1440    1333    1620    1920      10 

batch: 
   Length     Class      Mode 
       62 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
       62 character character 

}}
\keyword{datasets}
