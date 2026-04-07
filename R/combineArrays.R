#' Combine two RGChannelSet objects into one array object
#'
#' @param rgSet1 Input RGChannelSet object 1
#' @param rgSet2 Input RGChannelSet object 2
#' @param outType The array type of the combined RGChannelSet object
#'
#' @returns Output RGChannelSet object of the specified array type
#' @export
#'
#' @examples rgSetOutv1 <- combineArrays(rgSet1 = rgSetv1,
#' rgSet2 = rgSetv2) # input order does not matter here

combineArrays <- function(rgSet1, rgSet2,
                          outType = NULL){

    # Check to see if we have sensible inputs
    stopifnot(is(rgSet1, "RGChannelSet") & is(rgSet2, "RGChannelSet"))

    # Obtain annotations - if output type is not specified,
    anno1 <- BiocGenerics::annotation(rgSet1)[["array"]]
    anno2 <- BiocGenerics::annotation(rgSet2)[["array"]]
    options <- c("IlluminaHumanMethylation450k",
                 "IlluminaHumanMethylationEPIC",
                 "IlluminaHumanMethylationEPICv2")
    if(is.null(outType)) {
        outType <- options[options %in% c(anno1, anno2)][1]
    } else {
        outType <- match.arg(outType, options)
    }
    message(paste0("combining data as ", outType, "\n"))

    if(anno1 == anno2) {
        features <- intersect(rownames(rgSet1),
                              rownames(rgSet2))
        rgSet <- BiocGenerics::combine(rgSet1[features, ],
                                       rgSet2[features, ])
    } else if(!anno1 %in% options | !anno2 %in% options) {
        stop("These array types cannot be combined.")
    } else {
        if(anno1 != outType) {
            rgSet1 <- convertArray(rgSet1, outType = outType, verbose = verbose)
        }
        if(anno2 != outType) {
            rgSet2 <- convertArray(rgSet2, outType = outType, verbose = verbose)
        }
        features <- intersect(rownames(rgSet1),
                              rownames(rgSet2))
        rgSet <- BiocGenerics::combine(rgSet1[features, ],
                                       rgSet2[features, ])
        rgSet$ArrayTypes <- rep(
            x = c(anno1, anno2),
            times = c(ncol(rgSet1), ncol(rgSet2)))
        rgSet
    }
}
