# Combine reference and user's sample data sets together

Combine reference and user's sample data sets together

## Usage

``` r
combData(
  dataset,
  reference,
  class = c("rgset", "betas"),
  normType = c("Noob", "Funnorm", "Quantile", "Quantile.b", "None"),
  cellTypes = c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos",
    "Mono", "Neu", "NK", "Treg"),
  verbose = TRUE
)
```

## Arguments

- dataset:

  An RGchannelSet or a beta matrix of the users' DNA methylation data
  set

- reference:

  An RGchannelSet of the reference DNA methylation data set,
  preferentially the output from the
  [`getRef()`](https://maggie-fu.github.io/CellsPickMe/reference/getRef.md)
  function, to ensure the `pData()` formatting is appropriate

- class:

  A character specifying the argument `dataset`'s class, with the option
  or rgset (for RGChannelSet) or betas (for beta matrix)

- normType:

  A character specifying the method for normalizing user sample and
  reference data sets together, with the options of "Noob", "Funnorm",
  "Quantile", "Quantile.b", and "None" – refer to the Details section
  for more information on the options.

- cellTypes:

  A vector of characters specifying the cell types to be estimated. If
  the
  [`getRef()`](https://maggie-fu.github.io/CellsPickMe/reference/getRef.md)
  function was used, use the output as shown in the Example.

- verbose:

  A Boolean, specifying whether the function should be verbose

## Value

A list of dataframes, including the sample and reference data post
normalization, their pheno data, and the cell type information

## Examples

``` r
# Load example blood cell mixture, subsetted from the IDOL dataset (GSE110554)
test_dat <- CellsPickMe::IDOL_mixed_cells
#> Loading required namespace: minfi
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
# Obtain reference data set with the `getRef()` function
ref_dat <- getRef(ref = "IDOL", normType = "None")
#> 
#> see ?FlowSorted.Blood.EPIC and browseVignettes('FlowSorted.Blood.EPIC') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
# Combine sample and reference data sets together, followed by normalization (if selected)
comb_dat <- combData(dataset = test_dat, reference = ref_dat$reference, class = "rgset", normType = "None", cellTypes = ref_dat$cellTypes)
#> Combining Data with Flow Sorted Data and Normalizing.
#> Loading required package: IlluminaHumanMethylationEPICmanifest
#> Warning: there is no package called ‘IlluminaHumanMethylationEPICmanifest’
#> Error in getManifest(object): cannot load manifest package IlluminaHumanMethylationEPICmanifest
```
