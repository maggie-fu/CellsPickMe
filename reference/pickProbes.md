# Pick features for cell type prediction

Pick features for cell type prediction

## Usage

``` r
pickProbes(
  dataNormed,
  probeList = c("Ttest", "Caret", "IDOL", "DHS"),
  probeSelect = c("any", "both"),
  nProbes,
  caretMods = c("lasso", "EL", "BLR", "CART", "RF", "GBM", "GAnLDA", "GAnRF", "GAnNB",
    "GAnSVM", "GAnNN"),
  seed = 1,
  p.val = 0.05,
  min.delta.beta = 0,
  filterK = 1000,
  plotRef = TRUE,
  verbose = TRUE
)
```

## Arguments

- dataNormed:

  A list of dataframe containing the normalized data, output of
  [`combData()`](https://maggie-fu.github.io/CellsPickMe/reference/combData.md)

- probeList:

  A character, specifying how the probe should be selected. Options
  include "Ttest", "Caret_CV", "Caret_LOOCV", and "IDOL"

- probeSelect:

  If `probeSelect = Ttest`, specify how should the top probes be picked
  from T test output, options include "both", "any", "pval"

- nProbes:

  An integer specifying the number of probes to pick for each cell type

- caretMods:

  If `probeSelect %in% c("Caret_CV", "Caret_LOOCV")`, input a vector of
  characters that specify the models to use for feature selection

- seed:

  An integer specifying the seed for reproducibility

- p.val:

  If `probeSelect = Ttest`, specify a numeric value for maximum pvalue
  cutoff

- min.delta.beta:

  If `probeSelect = Ttest`, specify a numeric value for minimum delta
  beta (effect size) cutoff

- filterK:

  If `probeSelect %in% c("Caret_CV", "Caret_LOOCV")`, input an integer
  representing the number probes to input to the machine learning
  algorithms with top T test probes

- plotRef:

  A Boolean specifying whether to plot a heatmap showing the clustering
  performance using the probes selected

- verbose:

  A Boolean specifying whether the function should be verbose or not

## Value

A list containing the features selected by the specified methods, and
the coefficients of the selected probes for downstream prediction

## Examples
