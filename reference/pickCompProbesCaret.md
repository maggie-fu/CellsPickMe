# Pick features for cell type prediction with caret by repeated cross validation

Pick features for cell type prediction with caret by repeated cross
validation

## Usage

``` r
pickCompProbesCaret(
  betas,
  meta,
  ct,
  ps,
  min.delta.beta,
  p.val,
  caretMods,
  filterK = 1000,
  seed = 1234,
  plot = TRUE,
  verbose = TRUE
)
```

## Arguments

- betas:

  A beta matrix of reference DNA methylation data that will be used to
  select for features for cell type prediction

- meta:

  A data frame of phenotype data, specifying which of the reference
  samples are of what cell types

- ct:

  A vector of characters specifying the cell types to deconvolute

- ps:

  A character of either "any" or "both" to specify first line filter for
  T test prior to passing to machine learning algorithms

- min.delta.beta:

  A numeric variable defining T test minimum delta beta for first line
  filter prior to passing to machine learning algorithms

- p.val:

  A numeric variable defining maximum T test p.value for first line
  filter prior to passing to machine learning algorithms

- caretMods:

  A vector of characters, selecting the models to use to pick the cell
  type prediction features, options include "lasso", "EL", "BLR",
  "CART", "RF", "GBM", "PLDA", "GAnRF", "GAnNB", "GAnSVM", and "GAnNN"

- filterK:

  An integer, representing the number probes to input to the machine
  learning algorithms with top T test probes

- seed:

  An integer specifying the seed for reproducibility

- plot:

  A Boolean specifying whether to plot a heatmap showing the clustering
  performance using the probes selected

- verbose:

  A Boolean specifying whether the function should be verbose or not
