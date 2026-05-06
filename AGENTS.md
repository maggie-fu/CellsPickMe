## Cursor Cloud specific instructions

### Overview
CellsPickMe is an R package for predicting cell type proportions from DNA methylation (DNAme) data. It provides a pipeline from reference selection, to ML-based feature selection, to cell type prediction.

### Development commands
- **Load package in dev mode:** `Rscript -e 'devtools::load_all()'`
- **Install locally:** `Rscript -e 'devtools::install()'`
- **Generate docs:** `Rscript -e 'devtools::document()'`
- **Run R CMD check:** `Rscript -e 'devtools::check(args = "--no-manual")'`
- **Install deps:** `Rscript -e 'devtools::install_deps(dependencies = TRUE, upgrade = "never")'`

### Important notes
- There is **no test suite** (`tests/` directory does not exist). `devtools::test()` returns "No testing infrastructure found."
- `R CMD check` reports pre-existing errors/warnings in the package code itself (parse errors in examples, undeclared imports, LazyData compression). These are not environment issues.
- `getRef()` with UniBlood references downloads ~1GB files from Zenodo and is very slow. The IDOL reference (via ExperimentHub) is faster and gets cached locally.
- The package depends on Bioconductor packages. Install them via `BiocManager::install()`, not `install.packages()`.
- The `IlluminaHumanMethylation450kmanifest` and `IlluminaHumanMethylationEPICmanifest` packages are required at runtime but not listed in DESCRIPTION Imports/Suggests. They must be installed separately.
- System libraries required: `libcurl4-openssl-dev`, `libssl-dev`, `libxml2-dev`, `libfontconfig1-dev`, `libfreetype-dev`, `libpng-dev`, `libtiff-dev`, `libjpeg-dev`, `libharfbuzz-dev`, `libfribidi-dev`, `libgit2-dev`, `libuv1-dev`, `zlib1g-dev`, `libbz2-dev`, `liblzma-dev`, `libpcre2-dev`, `libgsl-dev`, `pandoc`.
