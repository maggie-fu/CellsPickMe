#' Convert an RGChannelSet into the outType array
#'
#' @param rgSet Input RGChannelSet object
#' @param outType The desired output array type, with the options of the 450k
#' array (IlluminaHumanMethylation450k), the EPICv1 array
#' (IlluminaHumanMethylationEPIC) or the EPICv2 array
#' (IlluminaHumanMethylationEPICv2)
#' @param verbose Verbosity of the function
#'
#' @returns Output RGChannelSet object converted to the specified array type
#' @export
#'
#' @examples rgSetv1 <- convertArray(rgSetv2,
#' outType = "IlluminaHumanMethylationEPIC",
#' verbose = T)

convertArray <- function(rgSet,
                         outType = c("IlluminaHumanMethylation450k",
                                     "IlluminaHumanMethylationEPIC",
                                     "IlluminaHumanMethylationEPICv2"),
                         verbose = verbose) {

    outType <- match.arg(outType)

    array <- BiocGenerics::annotation(rgSet)[["array"]]
    if(array == outType) stop("'rgSet' already in the 'outType' array type.")
    manifestIn <- .getManifest(array)
    manifestOut <- .getManifest(outType)
    annoOut <- if(outType == "IlluminaHumanMethylation450k"){
        "ilmn12.hg19"
    } else if(outType == "IlluminaHumanMethylationEPIC") {
        "ilm10b4.hg19"
    } else if(outType == "IlluminaHumanMethylationEPICv2") {
        "20a1.hg38"
    }

    keepAddresses <- list(
        I = NULL,
        II = NULL,
        SnpI = NULL,
        SnpII = NULL,
        Control = NULL)

    ### Probes of Type I
    probesIn <- getProbeInfo(manifestIn, type = "I")
    probesOut <- getProbeInfo(manifestOut, type = "I")
    if(array == "IlluminaHumanMethylationEPICv2") {
        # Remove duplicated probes
        probesIn$Name <- gsub("_.*", "", probesIn$Name)
        dupProbes <- table(probesIn$Name) %>%
            .[.>1] %>%
            names()
        probesIn <- probesIn[!probesIn$Name %in% dupProbes, ]
    }

    # Match the probes in the input and converted output
    commonNames <- intersect(probesIn$Name, probesOut$Name)
    probesIn <- probesIn[match(commonNames, probesIn$Name), ]
    probesOut <- probesOut[match(commonNames, probesOut$Name), ]

    # Remove probes with changed probe chemistry
    if(!all(probesIn$Color == probesOut$Color)){
        matchColor <- which(probesIn$Color == probesOut$Color)
        probesIn <- probesIn[matchColor, ]
        probesOut <- probesOut[matchColor, ]
    }

    # Sanity checks
    stopifnot(all(probesIn$Color == probesOut$Color))
    stopifnot(all(probesIn$ProbeSeqA == probesOut$ProbeSeqA))
    stopifnot(all(probesIn$ProbeSeqB == probesOut$ProbeSeqB))

    # Translating rgSetIn addresses to rgSetOut addresses
    translate <- c(probesOut$AddressA, probesOut$AddressB)
    names(translate) <- c(probesIn$AddressA, probesIn$AddressB)
    keepAddresses$I <- translate


    ### Probes of Type II
    probesIn <- getProbeInfo(manifestIn, type = "II")
    probesOut <- getProbeInfo(manifestOut, type = "II")
    if(array == "IlluminaHumanMethylationEPICv2") {
        # Remove duplicated probes
        probesIn$Name <- gsub("_.*", "", probesIn$Name)
        dupProbes <- table(probesIn$Name) %>%
            .[.>1] %>%
            names()
        probesIn <- probesIn[!probesIn$Name %in% dupProbes, ]
    }

    # Match the probes in the input and converted output
    commonNames <- intersect(probesIn$Name, probesOut$Name)
    probesIn <- probesIn[match(commonNames, probesIn$Name), ]
    probesOut <- probesOut[match(commonNames, probesOut$Name), ]

    # Remove probes with changed probe chemistry
    if(!all(probesIn$ProbeSeqA == probesOut$ProbeSeqA)){
        matchSeq <- which(probesIn$ProbeSeqA == probesOut$ProbeSeqA)
        probesIn <- probesIn[matchSeq, ]
        probesOut <- probesOut[matchSeq, ]
    }

    # Sanity checks
    stopifnot(all(probesIn$ProbeSeqA == probesOut$ProbeSeqA))

    # Translating rgSetIn addresses to rgSetOut addresses
    translate <- probesOut$AddressA
    names(translate) <- probesIn$AddressA
    keepAddresses$II <- translate


    ### Probes of Type SnpI
    probesIn <- getProbeInfo(manifestIn, type = "SnpI")
    probesOut <- getProbeInfo(manifestOut, type = "SnpI")
    if(array == "IlluminaHumanMethylationEPICv2") {
        # Remove duplicated probes
        probesIn$Name <- gsub("_.*", "", probesIn$Name)
        dupProbes <- table(probesIn$Name) %>%
            .[.>1] %>%
            names()
        probesIn <- probesIn[!probesIn$Name %in% dupProbes, ]
    }

    # Match the probes in the input and converted output
    commonNames <- intersect(probesIn$Name, probesOut$Name)
    probesIn <- probesIn[match(commonNames, probesIn$Name), ]
    probesOut <- probesOut[match(commonNames, probesOut$Name), ]

    # Sanity checks - sometimes the ProbeSeq A vs B is swapped - is that ok??
    stopifnot(all((probesIn$ProbeSeqA == probesOut$ProbeSeqB) |
                      (probesIn$ProbeSeqA == probesOut$ProbeSeqA)))
    stopifnot(all((probesIn$ProbeSeqB == probesOut$ProbeSeqB) |
                      (probesIn$ProbeSeqB == probesOut$ProbeSeqA)))

    # Translating rgSet2 addresses to rgSet1 addresses
    translate <- c(probesOut$AddressA, probesOut$AddressB)
    names(translate) <- c(probesIn$AddressA, probesIn$AddressB)
    keepAddresses$SnpI <- translate

    ### Probes of Type SnpII
    probesIn <- getProbeInfo(manifestIn, type = "SnpII")
    probesOut <- getProbeInfo(manifestOut, type = "SnpII")
    if(array == "IlluminaHumanMethylationEPICv2") {
        # Remove duplicated probes
        probesIn$Name <- gsub("_.*", "", probesIn$Name)
        dupProbes <- table(probesIn$Name) %>%
            .[.>1] %>%
            names()
        probesIn <- probesIn[!probesIn$Name %in% dupProbes, ]
    }

    # Match the probes in the input and converted output
    commonNames <- intersect(probesIn$Name, probesOut$Name)
    probesIn <- probesIn[match(commonNames, probesIn$Name),]
    probesOut <- probesOut[match(commonNames, probesOut$Name),]
    stopifnot(all(probesIn$ProbeSeqA == probesOut$ProbeSeqA))

    # Translating rgSet2 addresses to rgSet1 addresses
    translate <- probesOut$AddressA
    names(translate) <- probesIn$AddressA
    keepAddresses$SnpII <- translate

    ### Probes of Type Control
    probesIn <- getProbeInfo(manifestIn, type = "Control")
    probesOut <- getProbeInfo(manifestOut, type = "Control")
    commonAddress <- intersect(probesIn$Address, probesOut$Address)
    probesIn <- probesIn[match(commonAddress, probesIn$Address),]
    probesOut <- probesOut[match(commonAddress, probesOut$Address),]
    translate <- probesOut$Address
    names(translate) <- probesIn$Address
    keepAddresses$Control <- translate

    # Update rgSet

    keep <- unlist(unname(keepAddresses))
    rgSet <- rgSet[rownames(rgSet) %in% names(keep), ]
    rownames(rgSet) <- keep[rownames(rgSet)]
    annotation(rgSet) <- c(array = outType, annotation = annoOut)
    rgSet
}
