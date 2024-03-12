#' Calculate the stability of selected features with pvClust
#'
#' @import pvclust
#'
#' @param dataNormed A list of dataframe containing the normalized data, output of `combData()`
#' @param probes A list of probes to perform cell type prediction with and their corresponding coefficient, output of `pickProbes()`
#' @param parallel A SOCKcluster object specifying the number of processes to parallelize, created with `doParallel::makeCluster()`
#'
#' @return Approximately Unbiased (AU) measure of clusters based on hierarchical clustering results, representing the stability of each clusters via multiscale bootstrap resampling
#' @export
#'
#' @examples
#' # Load example blood cell mixture, subsetted from the IDOL dataset (GSE110554)
#' test_dat <- CellsPickMe::IDOL_mixed_cells
#' # Obtain reference data set with the `getRef()` function
#' ref_dat <- getRef(ref = "IDOL", normType = "None")
#' # Combine sample and reference data sets together, followed by normalization (if selected)
#' comb_dat <- combData(dataset = test_dat,
#' reference = ref_dat$reference, class = "rgset", normType = "None", cellTypes = ref_dat$cellTypes)
#' # Pick probes with repeated cross validation with T-test
#' probes <- pickProbes(dataNormed = comb_dat, probeList = "Ttest", probeSelect = "both", nProbes = 100, min.delta.beta = 0.05)
#' # Create parallelization clusters and calculate cluster stability
#' clustAU <- identClust(dataNormed = comb_dat, probes = probes, parallel = TRUE)


identClust <- function(dataNormed, probes , parallel = TRUE) {
    if(probes$probeList %in% c("Caret_CV", "Caret_LOOCV")) {
        out <- lapply(probes$coefs$probeCoefs, function(caretMod){
            clustout <- pvclust::pvclust(dataNormed$ref.n[rownames(caretMod), ], parallel = parallel)
            pvClust <- hc2split(clustout$hclust)$member
            ctClust <- sapply(dataNormed$cellType, function(x){
                return(which(dataNormed$refMeta$cellType == x))
            })
            clust <- lapply(ctClust, function(ct){
                y <- sapply(pvClust, function(x){identical(x, ct)}) %>% which(.)
                return(clustout$edges[y, ])
            })
            return(sapply(clust, function(x) x$au))
        })
    } else {
        clustout <- pvclust::pvclust(dataNormed$ref.n[rownames(probes$coefs), ], parallel = parallel)
        pvClust <- hc2split(clustout$hclust)$member
        ctClust <- sapply(comb_dat$cellType, function(x){
            return(which(comb_dat$refMeta$cellType == x))
        })
        clust <- lapply(ctClust, function(ct){
            y <- sapply(pvClust, function(x){identical(x, ct)}) %>% which(.)
            return(clustout$edges[y, ])
        })
        out <- sapply(clust, function(x) x$au)
    }
    return(out)
}
