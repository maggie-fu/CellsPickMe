
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
#> parallelly (1.37.0 -> 1.37.1) [CRAN]
#> digest     (0.6.34 -> 0.6.35) [CRAN]
#> data.table (1.15.0 -> 1.15.2) [CRAN]
#> package 'parallelly' successfully unpacked and MD5 sums checked
#> package 'digest' successfully unpacked and MD5 sums checked
#> package 'data.table' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\suanni\AppData\Local\Temp\RtmpMP1pqK\downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\suanni\AppData\Local\Temp\RtmpMP1pqK\remotes89f05b765c75\maggie-fu-CellsPickMe-e5a13a0/DESCRIPTION' ...  ✔  checking for file 'C:\Users\suanni\AppData\Local\Temp\RtmpMP1pqK\remotes89f05b765c75\maggie-fu-CellsPickMe-e5a13a0/DESCRIPTION'
#>       ─  preparing 'CellsPickMe': (518ms)
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>       ─  building 'CellsPickMe_0.0.0.9000.tar.gz'
#>      
#> 
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
#> Combining Data with Flow Sorted Data and Normalizing.
```

### Pick features that best distinguish cell type

The **CellsPickMe** package supports feature selection with either the
traditional T-test, or with machine-learning-based methods such as
elastic net and random forest to obtain a curated list of features that
are highly predictive of cell types

``` r
# Pick probes with T tests
probes <- pickProbes(dataNormed = comb_dat, 
                     probeList = "Ttest", #c("Ttest", "IDOL", "DHS")
                     probeSelect = "both", #c("both", "any", "pval")
                     nProbes = 100, # number of probes to pick for each cell type
                     p.val = 0.05,  # max pval
                     min.delta.beta = 0.05, # min delta beta
                     plotRef = T, # plot heatmap?
                     verbose = T)
#> Estimating Weights for Cell Type Prediction Based on Selected Probeset.
```

<img src="man/figures/README-pickProbes-1.png" width="100%" />

``` r

### Set up server for parallelization - run the code if picking probes with Caret
library(doParallel)
cl <- makeCluster(detectCores() - 1) # change as needed
registerDoParallel(cl)

# Pick probes with repeated cross validation with lasso and elastic net
probes <- pickProbes(dataNormed = comb_dat, 
                     probeList = "Caret_CV", #c("Caret_CV", "Caret_LOOCV")
                     caretMods = c("lasso", "EL"),  #c("lasso", "EL", "BLR", "CART", "RF", "GBM", "PLDA", "GAnRF", "GAnNB", "GAnSVM", "GAnNN")
                     filterK = 1000, # number of probes to put into the predictor for each cell type
                     seed = 1, 
                     plotRef = F, # plot heatmap?
                     verbose = F)
```

### Assess clustering stability

To further evaluate the performance of the selected probes, pvClust can
be applied to assess whether the picked probes can be used to generate
the correct cluster (cell type labeling) in reference data

``` r
clustAU <- identClust(dataNormed = comb_dat,
                      probes = probes,
                      parallel = TRUE)
#> Creating a temporary cluster...done:
#> socket cluster with 15 nodes on host 'localhost'
#> Multiscale bootstrap... Done.
#> Creating a temporary cluster...done:
#> socket cluster with 15 nodes on host 'localhost'
#> Multiscale bootstrap... Done.
#> Creating a temporary cluster...done:
#> socket cluster with 15 nodes on host 'localhost'
#> Multiscale bootstrap... Done.
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
                 verbose = F, 
                 cetygo = T) # CETYGO to assess reference appropriateness (RMSE evaluation)
```
