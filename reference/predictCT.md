# Estimate cell type proportion based on sample DNAm data and coefficients obtained from the reference data

Estimate cell type proportion based on sample DNAm data and coefficients
obtained from the reference data

## Usage

``` r
predictCT(
  dataNormed,
  probes,
  method,
  conditions = NULL,
  removenRBC = F,
  verbose = TRUE,
  cetygo = TRUE
)
```

## Arguments

- dataNormed:

  A list of dataframe containing the normalized data, output of
  [`combData()`](https://maggie-fu.github.io/CellsPickMe/reference/combData.md)

- probes:

  A list of probes to perform cell type prediction with and their
  corresponding coefficient, output of
  [`pickProbes()`](https://maggie-fu.github.io/CellsPickMe/reference/pickProbes.md)

- method:

  A character specifying the regression method. Options include "CP",
  "RPC", and "SVR"

- conditions:

  A character specifying specific conditions for model fitting. Default
  is NULL

- removenRBC:

  A Boolean specifying whether nucleated red blood cell (nRBC)
  proportion should be estimated, if using a reference with nRBC

- verbose:

  A Boolean specifying whether the function should be verbose or not

- cetygo:

  A Boolean specifying whether the CETYGO score should be calculated to
  estimate reference appropriateness

## Value

A matrix containing the estimated cell type proportion for the user
samples

## Examples

``` r
# Load example blood cell mixture, subsetted from the IDOL dataset (GSE110554)
test_dat <- CellsPickMe::IDOL_mixed_cells
# Obtain reference data set with the `getRef()` function
ref_dat <- getRef(ref = "IDOL", normType = "None")
#> see ?FlowSorted.Blood.EPIC and browseVignettes('FlowSorted.Blood.EPIC') for documentation
#> loading from cache
# Combine sample and reference data sets together, followed by normalization (if selected)
comb_dat <- combData(dataset = test_dat, reference = ref_dat$reference, class = "rgset", normType = "None", cellTypes = ref_dat$cellTypes)
#> Combining Data with Flow Sorted Data and Normalizing.
#> Loading required package: IlluminaHumanMethylationEPICmanifest
#> Warning: there is no package called ‘IlluminaHumanMethylationEPICmanifest’
#> Error in getManifest(object): cannot load manifest package IlluminaHumanMethylationEPICmanifest
# Pick probes with repeated cross validation with T-test
probes <- pickProbes(dataNormed = comb_dat, probeList = "Ttest", probeSelect = "both", nProbes = 100, min.delta.beta = 0.05)
#> Estimating Weights for Cell Type Prediction Based on Selected Probeset.
#> Error: object 'comb_dat' not found
# Estimate cell type proportion
out <- predictCT(dataNormed = comb_dat, probes = probes, method = "CP", cetygo = TRUE)
#> Estimating Composition Based on Selected Projection Method.
#> Error: object 'probes' not found
```
