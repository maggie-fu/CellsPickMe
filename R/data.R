#' IDOL mixed blood cell DNA methylation
#'
#' A subset of data from GSE110554
#'
#' @format ##
#' An RGChannelSet with 1,051,815 rows and 12 columns
#' This is a DNA methylation data set of artificially combined blood cell mixture, as measured on Illumina's EPIC microarray. The data represents the averaged DNA methylation level of the blood cell mixture.
#' The data is extracted from GSE110554, a dataset created and published by Salas et al. Please cite the following if using this data.
#' Salas LA, Koestler DC, Butler RA, Hansen HM et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome Biol 2018 May 29;19(1):64. PMID: 29843789
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554>
"IDOL_mixed_cells"


#' Middleton et al. Saliva reference pheno data
#'
#' The phenotype data of GSE147318
#'
#' @format ##
#' A dataframe of 60 rows and 1 column
#' This data contains the cell type label of the Middleton saliva reference, allowing for usability of the `getRef()` function when Middleton reference is called.
#' The data is extracted from GSE147318, a dataset created and published by Middleton et al. Please cite the following if using this data.
#' Middleton LYM, Dou J, Fisher J, Heiss JA et al. Saliva cell type DNA methylation reference panel for epidemiological studies in children. Epigenetics 2022 Jan-Feb;17(2):161-177. PMID: 33588693
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147318>
"Middleton_pd"


#' IDOL probes for 450k array - adult blood
#'
#' CpG names of probes selected by IDOL, for 450K DNAm microarray
#'
#' @format ##
#' A vector of characters of length 350
#' The data is obtained from a publication by Salas et al. Please cite the following if using this data.
#' Salas LA, Koestler DC, Butler RA, Hansen HM et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome Biol 2018 May 29;19(1):64. PMID: 29843789
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/29843789/>
"idol.a_450"

#' IDOL probes for EPIC array - adult blood
#'
#' CpG names of probes selected by IDOL, for EPIC DNAm microarray
#'
#' @format ##
#' A vector of characters of length 450
#' The data is obtained from a publication by Salas et al. Please cite the following if using this data.
#' Salas LA, Koestler DC, Butler RA, Hansen HM et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome Biol 2018 May 29;19(1):64. PMID: 29843789
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/29843789/>
"idol.a_EPIC"

#' IDOL probes for EPIC array - extended 12 cell type reference
#'
#' CpG names of probes selected by IDOL, for the extended blood cell type reference
#'
#' @format ##
#' A vector of characters of length 1200
#' The data is obtained from a publication by Salas et al. Please cite the following if using this data.
#' Salas LA, Zhang Z, Koestler DC, Butler RA et al. Enhanced cell deconvolution of peripheral blood using DNA methylation for high-resolution immune profiling. Nat Commun 2022 Feb 9;13(1):761. PMID: 35140201
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/35140201/>
"idol.a_ext"

#' IDOL probes for 450k array - cord blood
#'
#' CpG names of probes selected by IDOL, for cord blood cell type reference
#'
#' @format ##
#' A vector of characters of length 517
#' The data is obtained from a publication by Gervin et al.. Please cite the following if using this data.
#' Gervin, K., Salas, L. A., Bakulski, K. M., Van Zelm, M. C., Koestler, D. C., Wiencke, J. K., ... & Jones, M. J. (2019). Systematic evaluation and validation of reference and library selection methods for deconvolution of cord blood DNA methylation data. Clinical epigenetics, 11(1), 1-15.
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/29843789/>
"idol.c"
