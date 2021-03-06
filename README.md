# ProteoDisco

![license](https://img.shields.io/badge/license-GPL--3-blue.svg) [![GitHub issues](https://img.shields.io/github/issues/ErasmusMC-CCBC/ProteoDisco.svg)]() ![rversion](https://img.shields.io/badge/R%20version-%3E4.1.0-lightgrey.svg)

# Introduction

ProteoDisco is an R package to facilitate proteogenomics studies. 

It houses functions to create customized (variant) protein databases based on user-submitted genomic variants, splice-junctions, fusion genes and manual transcript sequences.
The flexible workflow can be adopted to suit a myriad of research and experimental settings.

## Citation

van de Geer, Wesley S. and van Riet, Job and van de Werken, Harmen J. G. [ProteoDisco: a flexible R approach to generate customized protein databases for extended search space of novel and variant proteins in proteogenomic studies](https://doi.org/10.1093/bioinformatics/btab809). *Bioinformatics* **2022** 38(5), 1437–1439.

# Installation

The latest development version can be installed directly from GitHub:

```R
# Require/install devtools package if not already installed.
if (!require("devtools")) install.packages("devtools", repos = "http://cran.r-project.org")

# Install ProteoDisco from GitHub.
devtools::install_github(repo = "ErasmusMC-CCBC/ProteoDisco")
library(ProteoDisco)
```

Or from BioConductor:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ProteoDisco")
```

# Usage

Please view the vignettes for instructions and tutorials on how to use this package.

```R
browseVignettes("ProteoDisco")
```
