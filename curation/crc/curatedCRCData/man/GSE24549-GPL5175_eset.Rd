\name{GSE24549-GPL5175_eset}
\alias{GSE24549-GPL5175_eset}
\docType{data}
\title{Transcriptome instability in colorectal cancer identified by exon microarray analyses: Associations with splicing factor expression levels and patient survival}
\description{Colorectal cancer is a heterogeneous disease molecularly characterized by inherent genomic instabilities, chromosome instability and microsatellite instability. In the present study we propose transcriptome instability as an analogue to genomic instability on the transcriptome level. Exon microarray data from two independent series of altoghether 160 colorectal cancer tissue samples was used for global alternative splicing detection using the FIRMA algorithm (aroma.affymetrix). The sample-wise amounts of these alternative splicing scores exceeding a defined threshold (deviating exon usage amounts) were summarized to provide the basis for description of transcriptome instability. This characteristic was shown to be associated with splicing factor expression levels and patient survival in both independent sample series.}
\usage{data(GSE24549-GPL5175_eset)}
\format{

  Authors: Sveen A, Agesen TH, Nesbakken A, Rognum TO et al.

  Lab: Oslo University Hospital, The Norwegian Radium Hospita.  Department of Cancer Prevention, Institute for Cancer Research.

  Contact: anita.sveen@rr-research.no

  Title: Transcriptome instability in colorectal cancer identified by exon microarray analyses: Associations with splicing factor expression levels and patient survival

  Url: http://pubget.com/search?q=21619627

  PubMedID: 21619627  
}
\details{
  Platform used: [HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]

  83 samples.

  GEO_platform_accession: GSE24549-GPL5175_eset
  platform_summary: GPL5175
  biomart_ID: [HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]

  Overall survival time-to-event summary (in years):
  Call: survfit(formula = Surv(time, cens) ~ -1)

  ___________________________________________________________________

  Available sample meta-data


  ___________________________________________________________________

  alt_sample_name:
  Length = 83.
  Class = character.
  Mode = character.

  sample_type:
  tumor = 83.

}
\source{
  http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE24549-GPL5175_eset
}
\keyword{datasets}