#' Combine reference and user's sample data sets together
#'
#' @param dataset An RGchannelSet or a beta matrix of the users' DNA methylation data set
#' @param reference An RGchannelSet of the reference DNA methylation data set, preferentially the output from the `getRef()` function, to ensure the `pData()` formatting is appropriate
#' @param class A character specifying the argument `dataset`'s class, with the option or rgset (for RGChannelSet) or betas (for beta matrix)
#' @param normType A character specifying the method for normalizing user sample and reference data sets together, with the options of "Noob", "Funnorm", "Quantile", "Quantile.b", and "None" -- refer to the Details section for more information on the options.
#' @param cellTypes A vector of characters specifying the cell types to be estimated. If the `getRef()` function was used, use the output as shown in the Example.
#' @param verbose A Boolean, specifying whether the function should be verbose
#'
#' @return A list of dataframes, including the sample and reference data post normalization, their pheno data, and the cell type information
#' @export
#'
#' @examples
#' # Load example blood cell mixture, subsetted from the IDOL dataset (GSE110554)
#' test_dat <- CellsPickMe::IDOL_mixed_cells
#' # Obtain reference data set with the `getRef()` function
#' ref_dat <- getRef(ref = "IDOL", normType = "None")
#' # Combine sample and reference data sets together, followed by normalization (if selected)
#' comb_dat <- combData(dataset = test_dat, reference = ref_dat$reference, class = "rgset", normType = "None", cellTypes = ref_dat$cellTypes)

combData <- function(dataset, reference, class = c("rgset", "betas"),
                     normType = c("Noob", "Funnorm", "Quantile", "Quantile.b", "None"),
                     cellTypes = c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem",
                                   "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg"),
                     verbose = TRUE) {

    if (!normType %in% c("Noob", "Funnorm", "Quantile", "Quantile.b", "None")) {
        stop("Please specify one of the available normalization methods")
    } else if (normType == "Funnorm") {
        processMethod <- "preprocessFunnorm"
        if (class == "betas") stop("A RGChannelSet is required for functional normalization")
    } else if (normType == "Noob") {
        processMethod <- "preprocessNoob"
        if (class == "betas") stop("A RGChannelSet or is required for Noob normalization")
    } else if (normType == "Quantile") {
        processMethod <- "preprocessQuantile"
        if (class == "betas") stop("A RGChannelSet is required for this quantile normalization method. If you only have betas, set normalization method to Quantile.b instead")
    } else if (normType == "Quantile.b") {
        if (class == "rgset") stop("Quantile.b is exclusively for beta matrix input. Use Quantile instead if you have an RGChannelSet")
    }

    sampCT <- rep("WBC", ncol(dataset))
    if (normType %in% c("Funnorm", "Noob", "Quantile")) {
        processMethod <- get(processMethod, envir = getNamespace('minfi'))
    }

    ### Data normalization
    if (verbose)
        cat("Combining Data with Flow Sorted Data and Normalizing.\n")
    reference <- reference[, minfi::pData(reference)$CellType %in% cellTypes]
    if (class == "betas") {
        # if (all(sampType %in% "CB")) {
        #     ref <- betas(reference[pData(reference)$SampType == "CB", ])
        # } else {
        #     ref <- betas(reference)
        # }
        if (class(reference) == "MethyLumiSet"){
            ref <- methylumi::betas(reference)
        } else if (class(reference) %in% c("MethylSet", "RGChannelSetExtended", "RGChannelSet", "GenomicRatioSet")){
            ref <- minfi::getBeta(reference)
        }
        commonprobe <- intersect(as.character(rownames(dataset)), as.character(rownames(ref)))
        if (normType == "Quantile.b") {
            comb <- cbind(dataset[commonprobe, ], ref[commonprobe, ])
            comb.n <- limma::normalizeQuantiles(comb)
            samp.n <- comb.n[, colnames(dataset)]
            ref.n <- comb.n[, colnames(ref)]
        } else if (normType == "None") {
            samp.n <- dataset[commonprobe, ]
            ref.n <- ref[commonprobe, ]
        }
    } else {
        dataset <- methods::as(dataset, "RGChannelSet")
        if (normType == "None") {
            if (class(reference) == "MethyLumiSet"){
                ref.n <- methylumi::betas(reference)
            } else if (class(reference) %in% c("RGChannelSetExtended", "RGChannelSet", "GenomicRatioSet")){
                ref.n <- minfi::getBeta(reference)
            }
            samp.n <- minfi::getBeta(dataset)
            commonprobe <- intersect(as.character(rownames(ref.n)), as.character(rownames(samp.n)))
            ref.n <- ref.n[commonprobe, ]
            samp.n <- samp.n[commonprobe, ]

        } else { # Combine all the datasets and normalize
            if(min(nrow(dataset), nrow(reference)) > 622400) {
                combRGset <- minfi::combineArrays(dataset, reference, outType = "IlluminaHumanMethylationEPIC")
            } else {
                combRGset <- minfi::combineArrays(dataset, reference, outType = "IlluminaHumanMethylation450k")
            }
            combRGset.N <- processMethod(combRGset)
            comb.n <- minfi::getBeta(combRGset.N)
            samp.n <- comb.n[, colnames(dataset)]
            ref.n <- comb.n[, colnames(reference)]
        }

    }
    combMeta <- data.frame(sampleNames = c(colnames(samp.n), colnames(ref.n)),
                           studyIndex = rep(c("user", "reference"), times = c(ncol(samp.n), ncol(ref.n))),
                           cellType = c(sampCT, as.character(minfi::pData(reference)$CellType)),
                           stringsAsFactors = FALSE)
    rownames(combMeta) <- combMeta$sampleNames
    refMeta <- combMeta[combMeta$studyIndex == "reference", ]
    sampMeta <- combMeta[combMeta$studyIndex == "user", ]
    return(list(samp.n = samp.n, ref.n = ref.n, refMeta = refMeta, sampMeta = sampMeta, cellTypes = cellTypes))
}
