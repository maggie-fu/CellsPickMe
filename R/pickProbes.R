#' Pick features for cell type prediction
#'
#' @importFrom rlang .data
#'
#' @param dataNormed A list of dataframe containing the normalized data, output of `combData()`
#' @param probeList A character, specifying how the probe should be selected. Options include "Ttest", "Caret_CV", "Caret_LOOCV", and "IDOL"
#' @param probeSelect If `probeSelect = Ttest`, specify how should the top probes be picked from T test output, options include "both", "any", "pval"
#' @param nProbes An integer specifying the number of probes to pick for each cell type
#' @param caretMods If `probeSelect %in% c("Caret_CV", "Caret_LOOCV")`, input a vector of characters that specify the models to use for feature selection
#' @param seed An integer specifying the seed for reproducibility
#' @param p.val If `probeSelect = Ttest`, specify a numeric value for maximum pvalue cutoff
#' @param min.delta.beta If `probeSelect = Ttest`, specify a numeric value for minimum delta beta (effect size) cutoff
#' @param filterK If `probeSelect %in% c("Caret_CV", "Caret_LOOCV")`, input an integer representing the number probes to input to the machine learning algorithms with top T test probes
#' @param plotRef A Boolean specifying whether to plot a heatmap showing the clustering performance using the probes selected
#' @param verbose A Boolean specifying whether the function should be verbose or not
#'
#' @return A list containing the features selected by the specified methods, and the coefficients of the selected probes for downstream prediction
#' @export
#'
#' @examples
#' # Load example blood cell mixture, subsetted from the IDOL dataset (GSE110554)
#' test_dat <- CellsPickMe::IDOL_mixed_cells
#' # Obtain reference data set with the `getRef()` function
#' ref_dat <- getRef(ref = "IDOL", normType = "None")
#' # Combine sample and reference data sets together, followed by normalization (if selected)
#' comb_dat <- combData(dataset = test_dat, reference = ref_dat$reference, class = "rgset", normType = "None", cellTypes = ref_dat$cellTypes)
#' # Pick probes with repeated cross validation with lasso and elastic net
#' probes <- pickProbes(dataNormed = comb_dat, probeList = "Caret_CV", caretMods = c("lasso", "EL"), filterK = 1000, seed = 1)


pickProbes <- function(dataNormed, probeList = c("Ttest", "Caret", "IDOL", "DHS"), probeSelect = c("any", "both"), nProbes,
                       caretMods = c("lasso", "EL", "BLR", "CART", "RF", "GBM", "GAnLDA", "GAnRF", "GAnNB", "GAnSVM", "GAnNN"),
                       seed = 1, p.val = 0.05, min.delta.beta = 0, filterK = 1000, plotRef = TRUE, verbose = TRUE) {
    ### Pick probes and estimate weights
    if (verbose)
        cat("Estimating Weights for Cell Type Prediction Based on Selected Probeset.\n")
    if (probeList == "Ttest") {
        coefs <- pickCompProbes2(betas = dataNormed$ref.n,
                                 meta = dataNormed$refMeta,
                                 ct = dataNormed$cellTypes,
                                 nP = nProbes,
                                 ps = probeSelect,
                                 trainingProbes = NULL,
                                 p.val = p.val,
                                 min.delta.beta = min.delta.beta,
                                 plot = plotRef) # Call the pickCompProbes2 function below to select the probes that can best discern cell types and calculate weights
    } else if (probeList == "Caret_CV") {
        coefs <- pickCompProbesCaret(betas = dataNormed$ref.n,
                                     meta = dataNormed$refMeta,
                                     ct = dataNormed$cellTypes,
                                     ps = probeSelect,
                                     p.val = p.val,
                                     min.delta.beta = min.delta.beta,
                                     caretMods = caretMods,
                                     verbose = verbose,
                                     filterK = filterK,
                                     plot = plotRef,
                                     seed = seed)
    } else if (probeList == "Caret_LOOCV") {
        coefs <- pickCompProbesCaretLOOCV(betas = dataNormed$ref.n,
                                          meta = dataNormed$refMeta,
                                          ct = dataNormed$cellTypes,
                                          ps = probeSelect,
                                          p.val = p.val,
                                          min.delta.beta = min.delta.beta,
                                          caretMods = caretMods,
                                          verbose = verbose,
                                          filterK = filterK,
                                          plot = plotRef,
                                          seed = seed)
    } else {
        if (probeList == "IDOL") {
            if ("nRBC" %in% dataNormed$cellTypes) {
                pLib <- idol.c
            } else if (nrow(dataNormed$samp.n) > 622399) {
                if ("CD8nv" %in% dataNormed$cellTypes) {
                    pLib <- idol.a_ext
                } else {
                    pLib <- idol.a_EPIC
                }
            } else {
                pLib <- idol.a_450
            }
        }
        # if (probeList == "DHS") {
        #     pLib <- DHS
        # }
        coefs <- pickCompProbes2(betas = dataNormed$ref.n,
                                 meta = dataNormed$refMeta,
                                 ct = dataNormed$cellTypes,
                                 nP = nProbes,
                                 trainingProbes = pLib,
                                 plot = plotRef)
    }

    return(list(probeList = probeList, coefs = coefs))
}
