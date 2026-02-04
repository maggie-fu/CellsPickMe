# Calculate the stability of selected features with pvClust

Calculate the stability of selected features with pvClust

## Usage

``` r
identClust(dataNormed, probes, parallel = TRUE)
```

## Arguments

- dataNormed:

  A list of dataframe containing the normalized data, output of
  [`combData()`](https://maggie-fu.github.io/CellsPickMe/reference/combData.md)

- probes:

  A list of probes to perform cell type prediction with and their
  corresponding coefficient, output of
  [`pickProbes()`](https://maggie-fu.github.io/CellsPickMe/reference/pickProbes.md)

- parallel:

  A SOCKcluster object specifying the number of processes to
  parallelize, created with `doParallel::makeCluster()`

## Value

Approximately Unbiased (AU) measure of clusters based on hierarchical
clustering results, representing the stability of each clusters via
multiscale bootstrap resampling

## Examples

``` r
# Load example blood cell mixture, subsetted from the IDOL dataset (GSE110554)
test_dat <- CellsPickMe::IDOL_mixed_cells
# Obtain reference data set with the `getRef()` function
ref_dat <- getRef(ref = "IDOL", normType = "None")
#> see ?FlowSorted.Blood.EPIC and browseVignettes('FlowSorted.Blood.EPIC') for documentation
#> loading from cache
# Combine sample and reference data sets together, followed by normalization (if selected)
comb_dat <- combData(dataset = test_dat,
reference = ref_dat$reference, class = "rgset", normType = "None", cellTypes = ref_dat$cellTypes)
#> Combining Data with Flow Sorted Data and Normalizing.
#> Loading required package: IlluminaHumanMethylationEPICmanifest
#> Warning: there is no package called ‘IlluminaHumanMethylationEPICmanifest’
#> Error in getManifest(object): cannot load manifest package IlluminaHumanMethylationEPICmanifest
# Pick probes with repeated cross validation with T-test
probes <- pickProbes(dataNormed = comb_dat, probeList = "Ttest", probeSelect = "both", nProbes = 100, min.delta.beta = 0.05)
#> Estimating Weights for Cell Type Prediction Based on Selected Probeset.
#> Error: object 'comb_dat' not found
# Create parallelization clusters and calculate cluster stability
clustAU <- identClust(dataNormed = comb_dat, probes = probes, parallel = TRUE)
#> Error: object 'probes' not found
```
