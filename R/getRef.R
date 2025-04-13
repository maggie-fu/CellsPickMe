#' Extract the correct reference data set
#'
#' @importFrom rlang .data
#'
#' @param ref A character representing the reference dataset, with the options of
#' "Reinius", "IDOL", "IDOL_extended", "Mixed", "Cord", "DLPFC", and "Middleton"
#' @param normType a character representing the normalization method, with the
#' options of "None", "Noob", "Funnorm", "Quantile"
#'
#' @return A list of the requested reference data set (RGchannelSet Object) and
#' the cell types of the sample (vector of character)
#' @export
#'
#' @examples
#' # Request the IDOL reference (2016) without normalization
#' getRef(ref = "IDOL", normType = "None")

getRef <- function(ref = c("Reinius", "IDOL", "Extended",
                           "UniBlood7", "UniBlood7",
                           "Cord", "DLPFC", "Middleton"),
                   normType = c("None", "Noob", "Funnorm", "Quantile")){

    if (!normType %in% c("Noob", "Funnorm", "Quantile", "None")) {
        stop("Please specify one of the available normalization methods")
    } else if (normType == "Funnorm") {
        processMethod <- "preprocessFunnorm"
    } else if (normType == "Noob") {
        processMethod <- "preprocessNoob"
    } else if (normType == "Quantile") {
        processMethod <- "preprocessQuantile"
    }

    if (normType != "None") {
        processMethod <- base::get(processMethod)
    }
    hub <- ExperimentHub::ExperimentHub()

    ### Reinus (adult blood) reference
    if(ref == "Reinius"){
        reference <- FlowSorted.Blood.450k::FlowSorted.Blood.450k
        cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
    }
    ### IDOL (adult blood) reference
    if(ref == "IDOL"){
        reference <- hub[["EH1136"]]
        minfi::pData(reference)$CellType <- minfi::pData(reference)$CellType %>%
            gsub("Neu", "Gran", .)
        cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
    }
    ### Extended (adult blood) reference
    if(ref == "Extended"){
        reference <- hub[["EH5425"]]
        n
        cellTypes <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv",
                       "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")
    }
    ### Cord (cord blood) reference
    if(ref == "Cord"){
        reference <- hub[["EH2256"]]
        minfi::pData(reference)$CellType <- minfi::pData(reference)$CellType %>%
            gsub("WholeBlood", "WB", .)
        cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")
    }
    ### UniBlood7 (adult+cord blood) reference
    if(ref == "UniBlood7"){
        message("Downloading UniBlood7.rda to temp dir and calling it.\n")
        dir <- tempdir()
        curl_download("https://zenodo.org/api/records/15204839/files/UniBlood7.rda/content",
                      paste0(dir, "/UniBlood7.rda"), quiet = FALSE)
        reference <- loadRData(paste0(dir, "/UniBlood7.rda"))
        cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")
    }
    ### UniBlood13 (adult+cord blood) reference
    if(ref == "UniBlood13"){
        message("Downloading UniBlood13.rda to temp dir and calling it.\n")
        dir <- tempdir()
        curl_download("https://zenodo.org/api/records/15204843/files/UniBlood13.rda/content",
                      paste0(dir, "/UniBlood13.rda"), quiet = FALSE)
        reference <- loadRData(paste0(dir, "/UniBlood13.rda"))
        cellTypes <- c("Bmem", "Bnv", "CD8mem", "CD8nv",
                       "CD4mem", "CD4nv", "Treg",
                       "Bas", "Eos", "Neu", "Mono", "NK", "nRBC")
    }
    ### UniBlood19 (adult+cord blood) reference
    if(ref == "UniBlood19"){
        message("Downloading UniBlood19.rda to temp dir and calling it.\n")
        dir <- tempdir()
        curl_download("https://zenodo.org/api/records/15204848/files/UniBlood19.rda/content",
                      paste0(dir, "/UniBlood19.rda"), quiet = FALSE)
        reference <- loadRData(paste0(dir, "/UniBlood19.rda"))
        cellTypes <- c("Bcell_cord", "Bmem", "Bnv",
                       "CD8T_cord", "CD8mem", "CD8nv",
                       "CD4T_cord", "CD4mem", "CD4nv", "Treg",
                       "Gran_cord", "Bas", "Eos", "Neu",
                       "Mono_cord", "Mono", "NK_cord", "NK", "nRBC")
    }
    ### DLPFC (brain) reference
    if(ref == "DLPFC"){
        reference <- FlowSorted.DLPFC.450k::FlowSorted.DLPFC.450k
        cellTypes <- c("NeuN_pos", "NeuN_neg")
    }
    ### Middleton (saliva) reference
    if(ref == "Saliva"){
        reference <- hub[["EH4539"]]
        minfi::pData(reference) <- Middleton_pd
        cellTypes <- c("epithelial", "immune")
    }

    ## Normalization
    if(normType != "None") {
        reference <- processMethod(reference)
    }

    return(list(reference = reference, cellTypes = cellTypes))
}
