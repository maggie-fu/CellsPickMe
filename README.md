
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CellsPickMe

<!-- badges: start -->
<!-- badges: end -->

**CellsPickMe** is a streamlined R package for cell type proportion
prediction in a heterogeneous tissue (i.e. deconvolution) based on DNA
methylation (DNAme) data. As cell identity is crucially linked to cell
function and DNAme profiles, estimation of cell type proportion is
instrumental to understanding major factors driving DNAme variability.
Accounting for inter-individual cellular heterogeneity is also vital in
epigenome-wide association studies to accurately evaluate the signals
associated with the variables of interest.

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
“Extended”, “UniBlood7”, “UniBlood13”, “UniBlood19”, “Cord”, “DLPFC”,
and “Middleton”), select one that is appropriate for your sample (based
on tissue and sample age) For more detail regarding the references,
check out the [Get Started](maggie-fu.github.io/CellsPickMe/) and
UniBlood Reference Generation page

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

Here is the T-test based feature selection method.

``` r
# Pick probes with T tests
probes <- pickProbes(dataNormed = comb_dat, 
                     probeList = "Ttest", #c("Ttest", "IDOL")
                     probeSelect = "both", #c("both", "any", "pval")
                     nProbes = 100, # number of probes to pick for each cell type
                     p.val = 0.05,  # max pval
                     min.delta.beta = 0.05, # min delta beta
                     plotRef = T, # plot heatmap?
                     verbose = T)
#> Estimating Weights for Cell Type Prediction Based on Selected Probeset.
```

<img src="man/figures/README-pickProbes1-1.png" width="100%" />

Alternatively, here are some options for machine-learning-based feature
selection.

``` r
### Set up server for parallelization - run the code if picking probes with Caret
library(doParallel)
cl <- makeCluster(detectCores() - 1) # change as needed
registerDoParallel(cl)

# Pick probes with repeated cross validation with lasso and elastic net (EN)
probes <- pickProbes(dataNormed = comb_dat, 
                     probeList = "Caret_CV", #c("Caret_CV", "Caret_LOOCV")
                     caretMods = c("RF", "EN"),  #c("lasso", "EN", "BLR", "CART", "RF", "GBM", "PLDA", "GAnRF", "GAnNB", "GAnSVM", "GAnNN")
                     probeSelect = "any",
                     p.val = 1, 
                     min.delta.beta = 0,
                     filterK = 1000, # number of probes to put into the predictor for each cell type
                     seed = 1, 
                     plotRef = F, # plot heatmap?
                     verbose = F)

# FilterK results used as input for feature selection
head(probes$coefs$probeList$CD8T$tTestTopK)
#> [1] "cg00219921" "cg25939861" "cg11531557" "cg18857618" "cg03388043"
#> [6] "cg02380585"
# Picked probes with RF and the estimated coefficients
head(probes$coefs$probeCoefs$RF)
#>                  CD8T      CD4T        NK      Bcell       Mono       Gran
#> cg00219921 0.13497027 0.8157856 0.9258168 0.95409930 0.94637918 0.94059532
#> cg25939861 0.08166750 0.8259936 0.6741324 0.91005096 0.90988840 0.91721765
#> cg11531557 0.10709924 0.5498983 0.9042382 0.93899463 0.95064388 0.94429583
#> cg02380585 0.92950484 0.9319455 0.0507419 0.02072647 0.02141933 0.02202468
#> cg02324835 0.07289746 0.3547489 0.8015559 0.91953938 0.92262876 0.91614711
#> cg27316811 0.15079427 0.4346930 0.9616850 0.95768196 0.97571470 0.97359512
# Picked probes with RF for each cell type and the cross validation performance 
head(probes$coefs$probeList$CD8T$RF)
#> $coefs
#>   [1] "cg00219921" "cg25939861" "cg11531557" "cg02380585" "cg02324835"
#>   [6] "cg27316811" "cg22512531" "cg07820882" "cg07661835" "cg04354393"
#>  [11] "cg10977115" "cg27111890" "cg26848126" "cg03357684" "cg21593001"
#>  [16] "cg22999502" "cg08777095" "cg07016730" "cg03591396" "cg24456340"
#>  [21] "cg04329870" "cg03172657" "cg15497761" "cg15880738" "cg12426092"
#>  [26] "cg23018755" "cg06329135" "cg27547947" "cg08212476" "cg04448700"
#>  [31] "cg18782488" "cg04759756" "cg14978069" "cg22193912" "cg04927910"
#>  [36] "cg12065750" "cg18222759" "cg01771871" "cg22807449" "cg00159523"
#>  [41] "cg23001918" "cg14374132" "cg00857031" "cg08563773" "cg22003402"
#>  [46] "cg25983560" "cg18174654" "cg07859138" "cg16563370" "cg25049210"
#>  [51] "cg01610979" "cg09234252" "cg11218434" "cg03998546" "cg09145126"
#>  [56] "cg25277116" "cg18425140" "cg14555759" "cg02360716" "cg08093026"
#>  [61] "cg04813880" "cg15925952" "cg27470208" "cg07836401" "cg21910220"
#>  [66] "cg07227647" "cg06131595" "cg07799277" "cg23257002" "cg18766861"
#>  [71] "cg19296671" "cg08704578" "cg15165154" "cg27619385" "cg26146569"
#>  [76] "cg05520031" "cg20048655" "cg15244101" "cg12118504" "cg12321040"
#>  [81] "cg09736953" "cg20867963" "cg12455066" "cg17671383" "cg19151402"
#>  [86] "cg26472802" "cg11067179" "cg01437221" "cg11796827" "cg07954091"
#>  [91] "cg01208318" "cg22208943" "cg01673307" "cg00188748" "cg23073974"
#>  [96] "cg04127242" "cg03899022" "cg01712812" "cg26779770" "cg14364151"
#> [101] "cg05872923" "cg10692325" "cg16932880" "cg02198513" "cg06943751"
#> [106] "cg17685731" "cg16338736" "cg07076915" "cg23893629" "cg27216684"
#> [111] "cg08087771" "cg13795831" "cg02470008" "cg16201233" "cg04281309"
#> [116] "cg16423339" "cg19482025" "cg03914662" "cg16235995" "cg14797308"
#> [121] "cg22264836" "cg11624270" "cg23314364" "cg26532208" "cg12916466"
#> [126] "cg16622355" "cg21912448" "cg09419354" "cg18814344" "cg05916684"
#> [131] "cg00482119" "cg16580362" "cg08940994" "cg24401262" "cg11093142"
#> [136] "cg01679017" "cg05249016" "cg16580955" "cg00083704" "cg09088496"
#> [141] "cg10680621" "cg18522931" "cg21340621" "cg16810279" "cg05280944"
#> [146] "cg16380497" "cg06766639" "cg03664994" "cg15979672" "cg23738183"
#> [151] "cg18349022" "cg08482894" "cg02061660" "cg17587388" "cg22481130"
#> [156] "cg12573049" "cg19235893" "cg17605068" "cg17729941" "cg16107927"
#> [161] "cg00133727" "cg21679873" "cg18653534" "cg07521825" "cg25195781"
#> [166] "cg08418872" "cg18300804" "cg15972391" "cg12616923" "cg15223598"
#> [171] "cg10114760" "cg02229211" "cg14140920" "cg09480603" "cg09978965"
#> [176] "cg05031116" "cg13988440" "cg16374333" "cg09588733" "cg10407106"
#> [181] "cg26157600" "cg02152068" "cg22581364" "cg02869300" "cg26363422"
#> [186] "cg04748191" "cg14585083" "cg01620165" "cg00714854" "cg12745332"
#> [191] "cg23594345" "cg23565155" "cg06689816" "cg18783296" "cg21720451"
#> [196] "cg22091063" "cg16315782" "cg10656896" "cg01867930" "cg11779949"
#> [201] "cg11788572" "cg12774429" "cg01305421" "cg24150662" "cg10254620"
#> [206] "cg06867417" "cg20168849" "cg22948808" "cg16681873" "cg04994519"
#> [211] "cg25504129" "cg01707663" "cg20282814" "cg09022607" "cg20239488"
#> [216] "cg01088404" "cg07466904" "cg26485825" "cg17061098" "cg03601289"
#> [221] "cg08769189" "cg23671196" "cg07674304" "cg11837532" "cg02888364"
#> [226] "cg04678315" "cg13564405" "cg11106946" "cg01592226" "cg17013990"
#> [231] "cg08326934" "cg02535635" "cg05976283" "cg04597393" "cg15345437"
#> [236] "cg18801945" "cg05496363" "cg12022722" "cg21267497" "cg03317021"
#> [241] "cg20988587" "cg12130067" "cg25535316" "cg06419846" "cg01002223"
#> [246] "cg25920106" "cg04212842" "cg26916319" "cg25096306" "cg26195363"
#> [251] "cg13974615" "cg08131309" "cg11537298" "cg03510285" "cg20367388"
#> [256] "cg06852314" "cg12183789" "cg12074856" "cg14468692" "cg21828461"
#> [261] "cg26055299" "cg03334072" "cg23487884" "cg17095147" "cg18854765"
#> [266] "cg06809693" "cg06494960" "cg23597825" "cg03456799" "cg00816039"
#> [271] "cg24182928" "cg26869829"
#> 
#> $performance
#>    Accuracy Kappa   Resample
#> 1         1     1 Fold4.Rep1
#> 2         1     1 Fold3.Rep2
#> 3         1     1 Fold2.Rep3
#> 4         1     1 Fold3.Rep1
#> 5         1     1 Fold2.Rep2
#> 6         1     1 Fold1.Rep3
#> 7         1     1 Fold2.Rep1
#> 8         1     1 Fold1.Rep2
#> 9         1     1 Fold5.Rep2
#> 10        1     1 Fold1.Rep1
#> 11        1     1 Fold5.Rep1
#> 12        1     1 Fold4.Rep2
#> 13        1     1 Fold3.Rep3
#> 14        1     1 Fold2.Rep4
#> 15        1     1 Fold1.Rep5
#> 16        1     1 Fold5.Rep5
#> 17        1     1 Fold1.Rep4
#> 18        1     1 Fold5.Rep4
#> 19        1     1 Fold4.Rep5
#> 20        1     1 Fold5.Rep3
#> 21        1     1 Fold4.Rep4
#> 22        1     1 Fold3.Rep5
#> 23        1     1 Fold4.Rep3
#> 24        1     1 Fold3.Rep4
#> 25        1     1 Fold2.Rep5
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

### Examine output

``` r
# Estimated proportions based on IDOL reference, Noob normalization, RF feature selection, and constraint projection
# CETYOGO score column shows the CETYGO score for each sample
head(out$RF)
#>                           CD8T       CD4T         NK      Bcell       Mono
#> 201868590193_R01C01 0.20697715 0.05755452 0.16832034 0.18737625 0.19889714
#> 201868590243_R02C01 0.05762550 0.15450029 0.03597759 0.05109316 0.07399635
#> 201868590267_R01C01 0.06765754 0.08656360 0.01636060 0.02983963 0.12079926
#> 201868590267_R05C01 0.12980382 0.15364182 0.03794995 0.02618161 0.09471185
#> 201869680008_R01C01 0.29098662 0.14490082 0.17136579 0.08515400 0.21821297
#> 201869680008_R03C01 0.09642377 0.13975148 0.04651692 0.03460012 0.08708938
#>                          Gran     CETYGO nCGmissing
#> 201868590193_R01C01 0.2138328 0.05366570          0
#> 201868590243_R02C01 0.6566618 0.03948246          0
#> 201868590267_R01C01 0.6989066 0.03682644          0
#> 201868590267_R05C01 0.5894684 0.04166677          0
#> 201869680008_R01C01 0.1258030 0.04975795          0
#> 201869680008_R03C01 0.6286077 0.04023896          0
head(out$EN)
#>                           CD8T       CD4T         NK      Bcell       Mono
#> 201868590193_R01C01 0.19518604 0.06283114 0.17031881 0.19263381 0.19788273
#> 201868590243_R02C01 0.04290187 0.15628262 0.04113237 0.05498541 0.07165128
#> 201868590267_R01C01 0.05427027 0.09287751 0.02320448 0.03369635 0.11590229
#> 201868590267_R05C01 0.11599535 0.15763322 0.04079190 0.02910951 0.09171784
#> 201869680008_R01C01 0.27794394 0.15197013 0.17364288 0.08785900 0.21793825
#> 201869680008_R03C01 0.08031931 0.14555998 0.04977482 0.03731340 0.08395700
#>                          Gran     CETYGO nCGmissing
#> 201868590193_R01C01 0.2114329 0.05232382          0
#> 201868590243_R02C01 0.6596637 0.03386079          0
#> 201868590267_R01C01 0.6994087 0.03209687          0
#> 201868590267_R05C01 0.5924770 0.03562926          0
#> 201869680008_R01C01 0.1220138 0.04378660          0
#> 201869680008_R03C01 0.6328493 0.03380049          0

# Basic visualization of deconvolution output
library(ggplot2)
library(dplyr)
library(reshape2)

plt <- out$RF[, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")] %>%
    as.data.frame(.) %>% 
    reshape2::melt()
ggplot(plt, aes(variable, value, color = variable)) + 
    geom_boxplot() + 
    geom_point(size = 2, position = position_jitter()) + 
    theme_bw() + 
    labs(x = "cell type", 
         y = "cell type proportion (%)", 
         title = "Blood cell type proportion") 
```

<img src="man/figures/README-visualization-1.png" width="100%" />

## Citation

The manuscript detailing CellsPickMe and its use is currently under
preparation. For more information about this please contact Maggie Fu at
<maggie.fu@bcchr.ca>.

## References

Depending on the options you used, please consider citing the following
references as this package is built on their data / code / papers.

1.  DS Vellame et al. (2023). Uncertainty quantification of
    reference-based cellular deconvolution algorithms. *Epigenetics* 18,
    1: 2137659. doi:
    [10.1080/15592294.2022.2137659](https://doi.org/10.1080/15592294.2022.2137659)

    - Please cite this paper if you set the parameter `cetygo = T` for
      `predictCT()`

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
