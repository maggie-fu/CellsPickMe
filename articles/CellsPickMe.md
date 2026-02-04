# CellsPickMe

### Introduction

**CellsPickMe** is a streamlined R package for cell type proportion
prediction in a heterogeneous tissue (i.e. deconvolution) based on DNA
methylation (DNAme) data. As cell identity is crucially linked to cell
function and DNAme profiles, estimation of cell type proportion is
instrumental to understanding major factors driving DNAme variability.
Accounting for inter-individual cellular heterogeneity is also vital in
epigenome-wide association studies to accurately evaluate the signals
associated with the variables of interest.

The **CellsPickMe** package takes DNAme data generated from Illumina
microarray and predict its cellular composition based on cell types
available in the reference profiles. *Notably, for each population of
**cells**, the package “**picks**” DNA**me** features that best predict
cellular identities with machine learning algorithms to improve
deconvolution performance.* Currently, the algorithm is compatible for
peripheral blood, cord blood, saliva, and brain (neuron vs non-neuron).
The table below illustrate the reference datasets available in the
**CellsPickMe** package for deconvolution. We curated the UniBlood
references (7, 13, and 19) to address the current challenge in
deconvoluting longitudinal and pediatric blood samples. Refer to the
\[UniBlood Reference Creation\]\[UniBlood Reference Creation\] article
for details of the UniBlood references.

##### Available reference datasets for tissues in CellsPickMe

[TABLE]

The figure below focused on our optimization of the deconvolution
pipeline in the blood tissue. The four main steps of cellular
deconvolution are reference selection, data normalization, feature
selection, and regression. The function that calls each of the step is
bolded on the left, and the detail of each function is described in the
[Usage](#usage) section.

CellsPickMe Pipeline, created with BioRender.com

## Cell Type Deconvolution

In addition to improved methods for feature selection to optimize
prediction process, the package also includes several features to assess
prediction accuracy.

The **CellsPickMe** package includes 5 core functions:

- [`getRef()`](https://maggie-fu.github.io/CellsPickMe/reference/getRef.md)
  finds and load the tissue-specific reference data set
- [`combData()`](https://maggie-fu.github.io/CellsPickMe/reference/combData.md)combines
  the reference and the user’s sample data sets together
- [`pickProbes()`](https://maggie-fu.github.io/CellsPickMe/reference/pickProbes.md)
  picks features that best distinguish cell types (multiple feature
  selection methods available)
- `clusterScore()` assesses the performance of selected features in
  clustering cell types in the reference
- [`predictCT()`](https://maggie-fu.github.io/CellsPickMe/reference/predictCT.md)
  estimate cell type proportions based on selected probes, with the
  ability to evaluate

The
[`pickProbes()`](https://maggie-fu.github.io/CellsPickMe/reference/pickProbes.md)
function can be computationally intensive, depending on the methods for
feature selection. Parallel computing is available and highly
recommended for improved efficiency. Detailed usage is available
[here](CellsPickMe::pickProbes())

The following data are required for performing the cell type proportion
prediction with **CellsPickMe**: - DNAme data (RGChannelSet Object is
recommended, beta matrix is accepted as well) \* An appropriate
reference can be obtained with the
[`getRef()`](https://maggie-fu.github.io/CellsPickMe/reference/getRef.md)
function. Otherwise the user can specify their own reference dataset as
well (needs to be a RGChannelSet Object with CellTypes information in
the pData)

------------------------------------------------------------------------

### Installation

``` r
library(devtools)
devtools::install_github("maggie-fu/CellsPickMe")
```

### Usage

The **CellsPickMe** package takes DNA methylation data generated from
Illumina microarray and predict its cellular composition based on cell
types available in the reference profiles. Currently, the algorithm is
compatible for peripheral blood, cord blood, saliva, and brain (neuron
vs non-neuron).

#### Obtain reference dataset

From the available list of reference datasets (“Reinius”, “IDOL”,
“Extended”, “Cord”, “UniBlood7”, “UniBlood13”, “UniBlood19”, “DLPFC”,
and “Middleton”), select one that is appropriate for your sample (based
on tissue and sample age)

``` r
library(CellsPickMe)

# Request the UniBlood19 reference with no normalization
ref_dat <- getRef(ref = "IDOL", normType = "None")
```

#### Normalize sample and reference datasets together

Normalize user’s sample and reference data sets together to reduce batch
effect and improve prediction accuracy. Option for `normType` includes
“Noob”, “Funnorm”, “Quantile”, “Quantile.b”, “None”, with the first 3
options being exclusively for RGChannelSet objects, and Quantile.b for
beta matrix.

``` r
# Load example blood cell mixture, subsetted from the IDOL dataset (GSE110554)
test_dat <- CellsPickMe::IDOL_mixed_cells
 
# Combine sample and reference data sets together, followed by normalization (if selected)
comb_dat <- combData(dataset = test_dat, 
                     reference = ref_dat$reference, 
                     cellTypes = ref_dat$cellTypes, 
                     class = "rgset",               #c("rgset", "betas")
                     normType = "None")             # c("Noob", "Funnorm", "Quantile", "Quantile.b", "None")
```

#### Pick features that best distinguish cell type

The **CellsPickMe** package supports feature selection with either the
traditional T-test, or with machine-learning-based methods such as
elastic net and random forest to obtain a curated list of features that
are highly predictive of cell types.

[TABLE]

``` r
# Pick probes with T tests
probes <- pickProbes(dataNormed = comb_dat, 
                     probeList = "Ttest", #c("Ttest", "Caret_CV", "Caret_LOOCV", "IDOL")
                     probeSelect = "both", #c("both", "any")
                     nProbes = 100, # number of probes to pick for each cell type
                     p.val = 0.05,  # max pval
                     min.delta.beta = 0.05, # min delta beta
                     plotRef = T, # plot heatmap?
                     verbose = T)

### Set up server for parallelization - run the code if picking probes with Caret
library(doParallel)
cl <- makeCluster(detectCores() - 1) # change as needed
registerDoParallel(cl)

# Pick probes with repeated cross validation with lasso and elastic net
probes <- pickProbes(dataNormed = comb_dat, 
                     probeList = "Caret_LOOCV", #c("Ttest", "Caret_CV", "Caret_LOOCV", "IDOL")
                     caretMods = c("lasso", "EN"),  #c("lasso", "EN", "BLR", "CART", "RF", "GBM")
                     filterK = 1000, # number of probes to put into the predictor for each cell type
                     seed = 1, 
                     plotRef = F, # plot heatmap?
                     verbose = F)

# FilterK results used as input for feature selection
head(probes$coefs$probeList$CD8T$tTestTopK)

# Picked probes with EN and the estimated coefficients
head(probes$coefs$probeCoefs$EN)

# Picked probes with EN for each cell type and the cross validation performance 
head(probes$coefs$probeList$CD8mem$EN$coefs)
```

#### Assess clustering stability

To further evaluate the performance of the selected probes, pvClust can
be applied to assess whether the picked probes can be used to generate
the correct cluster (cell type labeling) in reference data

``` r
clustAU <- identClust(dataNormed = comb_dat,
                      probes = probes,
                      parallel = TRUE)
```

#### Estimate cell type proportion

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

### Citation

The manuscript detailing CellsPickMe and its use is currently under
preparation. For more information about this please contact Maggie Fu at
<maggie.fu@bcchr.ca>.

### References

Depending on the options you used, please consider citing the following
references as this package is built on their data / code / papers.

1.  DS Vellame et al. (2023). Uncertainty quantification of
    reference-based cellular deconvolution algorithms. *Epigenetics* 18,
    1: 2137659. doi:
    [10.1080/15592294.2022.2137659](https://doi.org/10.1080/15592294.2022.2137659)

    - Please cite this paper if you set the parameter `cetygo = T` for
      [`predictCT()`](https://maggie-fu.github.io/CellsPickMe/reference/predictCT.md)

2.  LE Reinius et al. (2012). Differential DNA methylation in purified
    human blood cells: implications for cell lineage and studies on
    disease susceptibility. *PloS one*. *7*(7), e41361. doi:
    [10.1371/journal.pone.0041361](https://doi.org/10.1371/journal.pone.0041361)

    - Please cite this paper if you used this reference dataset,
      i.e. `getRef(ref = "Reinius")`

3.  DC Koestler et al. (2016). Improving cell mixture deconvolution by
    identifying optimal DNA methylation libraries (IDOL). *BMC
    bioinformatics*. 17, 120. doi:
    [10.1186/s12859-016-0943-7](https://doi.org/10.1186/s12859-016-0943-7)

    - Please cite this paper if you used this reference dataset,
      i.e. `getRef(ref = "IDOL")`

4.  LA Salas et al. (2022) Enhanced cell deconvolution of peripheral
    blood using DNA methylation for high-resolution immune profiling.
    *Nat Commun*. 13, 761. doi:
    [10.1038/s41467-021-27864-7](https://doi.org/10.1038/s41467-021-27864-7)

    - Please cite this paper if you used this reference dataset,
      i.e. `getRef(ref = "Extended")`

5.  K Gervin et al. (2019). Systematic evaluation and validation of
    reference and library selection methods for deconvolution of cord
    blood DNA methylation data. *Clinical epigenetics*. 11, 1-15. doi:
    [10.1186/s13148-019-0717-y](https://doi.org/10.1186/s13148-019-0717-y)

    - Please cite this paper if you used this reference dataset,
      i.e. `getRef(ref = "Cord")`

6.  LY Middleton et al. (2022). Saliva cell type DNA methylation
    reference panel for epidemiological studies in children.
    *Epigenetics*, *17*(2), 161-177. doi:
    [10.1371/journal.pone.0041361](https://doi.org/10.1371/journal.pone.0041361)

    - Please cite this paper if you used this reference dataset,
      i.e. `getRef(ref = "Middleton")`

7.  J Guintivano et al. (2013). A cell epigenotype specific model for
    the correction of brain cellular heterogeneity bias and its
    application to age, brain region and major depression.
    *Epigenetics*, *8*(3), 290-302. doi:
    [10.4161/epi.23924](https://doi.org/10.4161/epi.23924)

    - Please cite this paper if you used this reference dataset,
      i.e. `getRef(ref = "DLPFC")`

8.  TJ Triche, et al. (2013). Low-level processing of Illumina Infinium
    DNA Methylation BeadArrays. *Nucleic Acids Res*. 41, e90. doi:
    [10.1093/nar/gkt090](http://www.dx.doi.org/10.1093/nar/gkt090).

    - Please cite this paper if you used this normalization method,
      i.e. `getRef(normType = "Noob")` or `combData(normType = "Noob")`

9.  JP Fortin et al. (2014). Functional normalization of 450k
    methylation array data improves replication in large cancer studies.
    *Genome Biology* 15, 503. doi:
    [10.1186/s13059-014-0503-2](http://www.dx.doi.org/10.1186/s13059-014-0503-2).

    - Please cite this paper if you used this normalization method,
      i.e. `getRef(normType = "Funnorm")` or
      `combData(normType = "Funnorm")`

10. N Touleimat and J Tost. (2012). *Complete pipeline for Infinium
    Human Methylation 450K BeadChip data processing using subset
    quantile normalization for accurate DNA methylation estimation.*
    *Epigenomics* 4, 325-341. doi:
    [10.2217/epi.12.21](https://doi.org/10.2217/epi.12.21)

    - Please cite this paper if you used this normalization method,
      i.e. `getRef(normType = "Quantile")` or
      `combData(normType = "Quantile")`
