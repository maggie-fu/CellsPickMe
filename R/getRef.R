#' Extract the correct reference data set
#'
#' @param ref A character representing the reference dataset, with the options of "Reinius", "IDOL", "IDOL_extended", "Mixed", "Cord", "DLPFC", and "Middleton"
#' @param normType a character representing the normalization method, with the options of "None", "Noob", "Funnorm", "Quantile"
#'
#' @return A list of the requested reference data set (RGchannelSet Object) and the cell types of the sample (vector of character)
#' @export
#'
#' @examples
#' # Request the IDOL extended reference (2022) without normalization
#' getRef(ref = "IDOL_extended", normType = "None")
#' # Request the Middleton reference (2022) normalized with functional normalization
#' getRef(ref = "Middleton", normType = "Funnorm")
#'
getRef <- function(ref = c("Reinius", "IDOL", "IDOL_extended", "Mixed",
                           "Cord", "DLPFC", "Middleton"),
                   normType = c("None", "Noob", "Funnorm", "Quantile")){

    hub <- ExperimentHub::ExperimentHub()

    ### Reinus (adult blood) reference
    if(ref == "Reinius"){
        library(FlowSorted.Blood.450k)
        reference <- utils::data(FlowSorted.Blood.450k)
        cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
    }
    ### IDOL (adult blood) reference
    if(ref == "IDOL"){
        reference <- hub[["EH1136"]]
        Biobase::pData(reference)$CellType <- Biobase::pData(reference)$CellType %>%
            gsub("Neu", "Gran", .)
        cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
    }
    ### Extended (adult blood) reference
    if(ref == "IDOL_extended"){
        reference <- hub[["EH5425"]]
        cellTypes <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv",
                       "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")
    }
    ### Cord (cord blood) reference
    if(ref == "Cord"){
        reference <- hub[["EH2256"]]
        Biobase::pData(reference)$CellType <- Biobase::pData(reference)$CellType %>%
            gsub("WholeBlood", "WB", .)
        cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")
    }
    ### Mixed (adult+cord blood) reference
    if(ref == "Mixed"){
        reference <- base::readRDS("data/Mixed_Reference_RGset.rds")
        cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")
    }
    ### DLPFC (brain) reference
    if(ref == "DLPFC"){
        library(FlowSorted.DLPFC.450k)
        reference <- FlowSorted.DLPFC.450k
        cellTypes <- c("NeuN_pos", "NeuN_neg")
    }
    ### Middleton (saliva) reference
    if(ref == "Saliva"){
        reference <- hub[["EH4539"]]
        Middleton_pd <- utils::read.table("data/Middleton_pd.txt")
        Biobase::pData(reference) <- Middleton_pd
        cellTypes <- c("epithelial", "immune")
    }
    return(list(reference = reference, cellTypes = cellTypes))
}
