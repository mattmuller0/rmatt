# rmatt

Install the package using the following command:

```R
devtools::install_github("mattmuller0/rmatt")
```

## General Functions

- general_functions.R: General functions
- plotting_functions.R: Functions for plotting (mainly ggplot2)
- utils.R: Utility functions
- exports.R: Exports of some common functions from other packages

## Omics Functions

- converting_functions.R: Functions for converting data and annotations
- filtering_functions.R: Functions for filtering rna-seq data
- rna_functions.R: Functions for RNA analysis
- proteomics_functions.R: Functions for proteomics analysis (targeted proteomics mainly)
- phenomics_functions.R: Functions for biobank phenomics analysis
- enrichment_functions.R: Functions for enrichment analysis
- signature_functions.R: Functions for -omics signature analysis

## Stats Functions

- stats_functions.R: Functions for statistics and survival analysis
- clustering_functions.R: Functions for clustering analysis
- survival_functions.R: Functions for survival analysis

## Installing Dependencies

This package depends on several CRAN and Bioconductor packages. To install all required dependencies, run the following in R:

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# List of required packages
cran_packages <- c(
  "cowplot", "dplyr", "ggplot2", "ggpubr", "glue", "purrr", "stringr", "tibble", "tidyr", "tidyverse", "broom", "forcats", "ggbiplot", "ggrepel", "ggsurvfit", "glmnet", "rlang", "rstatix", "survival", "magrittr", "RColorBrewer", "enrichR", "tableone"
)

bioc_packages <- c(
  "AnnotationDbi", "BiocManager", "ComplexHeatmap", "DESeq2", "NMF", "OlinkAnalyze", "S4Vectors", "SummarizedExperiment", "circlize", "clusterProfiler", "edgeR", "enrichplot", "ggtree", "limma", "msigdbr", "org.Hs.eg.db", "singscore"
)

# Install CRAN packages
install.packages(setdiff(cran_packages, rownames(installed.packages())))

# Install Bioconductor packages
BiocManager::install(setdiff(bioc_packages, rownames(installed.packages())))
```
