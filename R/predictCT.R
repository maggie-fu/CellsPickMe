#' Estimate cell type proportion based on sample DNAm data and coefficients obtained from the reference data
#'
#' @param dataNormed A list of dataframe containing the normalized data, output of `combData()`
#' @param probes A list of probes to perform cell type prediction with and their corresponding coefficient, output of `pickProbes()`
#' @param method A character specifying the regression method. Options include "CP", "RPC", and "SVR"
#' @param conditions A character specifying specific conditions for model fitting. Default is NULL
#' @param removenRBC A Boolean specifying whether nucleated red blood cell (nRBC) proportion should be estimated, if using a reference with nRBC
#' @param verbose A Boolean specifying whether the function should be verbose or not
#' @param cetygo A Boolean specifying whether the CETYGO score should be calculated to estimate reference appropriateness
#'
#' @return A matrix containing the estimated cell type proportion for the user samples
#' @export
#'
#' @examples
#' # Load example blood cell mixture, subsetted from the IDOL dataset (GSE110554)
#' test_dat <- CellsPickMe::IDOL_mixed_cells
#' # Obtain reference data set with the `getRef()` function
#' ref_dat <- getRef(ref = "IDOL", normType = "None")
#' # Combine sample and reference data sets together, followed by normalization (if selected)
#' comb_dat <- combData(dataset = test_dat, reference = ref_dat$reference, class = "rgset", normType = "None", cellTypes = ref_dat$cellTypes)
#' # Pick probes with repeated cross validation with T-test
#' probes <- pickProbes(dataNormed = comb_dat, probeList = "Ttest", probeSelect = "both", nProbes = 100, min.delta.beta = 0.05)
#' # Estimate cell type proportion
#' out <- predictCT(dataNormed = comb_dat, probes = probes, method = "CP", cetygo = TRUE)

predictCT <- function(dataNormed, probes, method, conditions = NULL, removenRBC = F, verbose = TRUE, cetygo = TRUE){
    if (verbose)
        cat("Estimating Composition Based on Selected Projection Method.\n")
    projectionMethod <- get(method)

    if (probes$probeList %in% c("Caret_LOOCV", "Caret_CV")){
        out <- lapply(probes$coefs$probeCoefs, function(coefs){
            if(nrow(coefs) <= ncol(coefs)) {
                message("There are fewer features than the number of cell types you're trying to predict! skipping this prediction")
            } else {
                mat <- dataNormed$samp.n[rownames(coefs), ]
                mat <- as.matrix(mat[stats::complete.cases(mat), ])
                coefs <- as.matrix(coefs)
                if(removenRBC){
                    coefs <- coefs[, colnames(coefs) != "nRBC"] # For peripheral blood, remove nRBC in prediction
                }
                if (ncol(coefs) > 0) {
                    if (is.null(conditions)) counts <- projectionMethod(samp.n = mat, coef = coefs) %>%
                            as.data.frame() # Using the weights generated in the last step to stimate the proportion of each cell type
                    else counts <- projectionMethod(samp.n = mat, coef = coefs, conditions = conditions) %>%
                            as.data.frame()
                }
                # counts[counts < 0] <- 0
                out <- counts

                if(cetygo){
                    YIN <- dataNormed$samp.n[rownames(coefs), ]
                    CETYGO <- sapply(seq_len(nrow(counts)), function(x) {
                        getErrorPerSample(applyIndex = x,
                                          predictedIN = counts,
                                          coefDataIN = coefs,
                                          betasBulkIN = YIN)
                    })
                    nMissingAll <- nrow(coefs) - nrow(YIN)
                    nCGmissing <- apply(YIN, 2, function(x) {
                        sum(is.na(x)) + nMissingAll
                    })
                    out <- cbind(counts, CETYGO, nCGmissing)
                }
                return(out)
            }
        })
    } else {
        mat <- dataNormed$samp.n[rownames(probes$coefs), ]
        mat <- as.matrix(mat[stats::complete.cases(mat), ])
        coefs <- as.matrix(probes$coefs)
        if(removenRBC){
            coefs <- coefs[, colnames(coefs) != "nRBC"] # For peripheral blood, remove nRBC in prediction
        }
        if (ncol(mat) > 0) {
            if (is.null(conditions)) counts <- projectionMethod(samp.n = mat, coef = coefs) %>%
                    as.data.frame()
            else counts <- projectionMethod(samp.n = mat, coef = coefs, conditions = conditions) %>%
                    as.data.frame()
        }
        out <- counts

        if(cetygo){
            YIN <- dataNormed$samp.n[rownames(coefs), ]
            CETYGO <- sapply(seq_len(nrow(counts)), function(x) {
                getErrorPerSample(applyIndex = x,
                                  predictedIN = counts,
                                  coefDataIN = coefs,
                                  betasBulkIN = YIN)
            })
            nMissingAll <- nrow(probes$coefs) - nrow(YIN)
            # nCGmissing <- apply(YIN, 2, function(x) {
            #     sum(is.na(x)) + nMissingAll
            # })
            out <- cbind(counts, CETYGO)
        }
    }
    return(out)
}


