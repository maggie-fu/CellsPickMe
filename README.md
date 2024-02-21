
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CellsPickMe

<!-- badges: start -->
<!-- badges: end -->

**CellsPickMe** is a comprehensive R package for cell type proportion
prediction based on DNA methylation data. As cell identity is crucially
linked to cell function and DNA methylation profile, estimation of cell
type proportion is instrumental to understanding major factors driving
DNA methylation variability. Accounting for inter-individual cellular
heterogeneity is also vital in epigenome-wide association studies (EWAS)
to accurately evaluate the signals associated with the variables of
interest.

## Installation

You can install the development version of CellsPickMe from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
devtools::install_github("maggie-fu/CellsPickMe")
```

## Usage

The **CellsPickMe** package takes DNA methylation data generated from
Illumina microarray and predict its cellular composition based on cell
types available in the reference profiles. Currently, the algorithm is
compatible for peripheral blood, cord blood, saliva, and brain (neuron
vs non-neuron).

### Obtain reference dataset

From the available list of reference datasets (“Reinius”, “IDOL”,
“IDOL_extended”, “Mixed”, “Cord”, “DLPFC”, and “Middleton”), select one
that is appropriate for your sample (based on tissue and sample age)

``` r
library(CellsPickMe)

# Request the IDOL reference (2016) with no normalization
ref_dat <- getRef(ref = "IDOL", normType = "None")
```

### Normalize sample and reference datasets together

Normalize user’s sample and reference data sets together to reduce batch
effect and improve prediction accuracy

``` r
# Load example blood cell mixture, subsetted from the IDOL dataset (GSE110554)
test_dat <- CellsPickMe::IDOL_mixed_cells
 
# Combine sample and reference data sets together, followed by normalization (if selected)
comb_dat <- combData(dataset = test_dat, 
                     reference = ref_dat$reference, 
                     class = "rgset", 
                     normType = "Noob", 
                     cellTypes = ref_dat$cellTypes)
```

### Pick features that best distinguish cell type

The **CellsPickMe** package supports feature selection with either the
traditional T-test, or with machine-learning-based methods such as
elastic net and random forest to obtain a curated list of features that
are highly predictive of cell types

``` r
# Pick probes with T tests
probes <- pickProbes(dataNormed = comb_dat, 
                     probeList = "Ttest", #c("Ttest", "Caret", "IDOL", "DHS")
                     probeSelect = "both", #c("both", "any", "pval")
                     nProbes = 100, # number of probes to pick for each cell type
                     p.val = 0.05,  # max pval
                     min.delta.beta = 0.05, # min delta beta
                     plotRef = F, # plot heatmap?
                     verbose = T)

### Set up server for parallelization - run the code if picking probes with Caret or doing pvClust
library(doParallel)
cl <- makeCluster(10) # change as needed
registerDoParallel(cl)

# Pick probes with repeated cross validation with lasso and elastic net
probes <- pickProbes(dataNormed = comb_dat, 
                     probeList = "Caret_CV", #c("Ttest", "Caret_CV", "Caret_LOOCV", "IDOL", "DHS")
                     caretMods = c("lasso", "EL"),  #c("lasso", "EL", "BLR", "CART", "RF", "GBM", "PLDA", "GAnRF", "GAnNB", "GAnSVM", "GAnNN")
                     filterK = 1000, # number of probes to put into the predictor for each cell type
                     seed = 1, 
                     plotRef = T, # plot heatmap?
                     verbose = T)
```

### Assess clustering stability

To further evaluate the performance of the selected probes, pvClust can
be applied to assess whether the picked probes can be used to generate
the correct cluster (cell type labeling) in reference data

``` r
clustout <- pvclust(comb_dat$ref.n[rownames(probes$coefs$probeCoefs$lasso), ], parallel = cl)
pvClust <- hc2split(clustout$hclust)$member
ctClust <- sapply(ref$cellTypes, function(x){ # can be replaced with user cellTypes arguments
    return(which(comb_dat$refMeta$cellType == x))
}) 
identClust <- lapply(ctClust, function(ct){
    y <- sapply(pvClust, function(x){identical(x, ct)}) %>% which(.)
    return(clustout$edges[y, ])
})
```

### Estimate cell type proportion

Finally, the selected features are used to estimate cell type
proportions in the sample data set. We also incorporate the CETYGO score
(see Reference) to estimate prediction performance even in the absence
of a validation cohort / ground truth cell count.

``` r
out <- predictCT(dataNormed = comb_dat, 
                 probes = probes, 
                 method = "CP",  #c("CP", "RPC", "SVR")
                 removenRBC = F, # remove nRBC?
                 verbose = T, 
                 cetygo = T) # CETYGO to assess reference appropriateness (RMSE evaluation)
```
