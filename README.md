# scMultiMap

`scMultiMap` is an R package for inferring cell-type-specific peak-gene associations using single-cell multimodal data. It implements the statistical method proposed in article [scMultiMap: Cell-type-specific mapping of enhancers and target genes from single-cell multimodal data](https://doi.org/10.1038/s41467-025-59306-z).

**Citation**: Chang Su, Dongsoo Lee, Peng Jin and Jingfei Zhang. (2025). scMultiMap: Cell-type-specific mapping of enhancers and target genes from single-cell multimodal data. *Nature Communications*.


## Installation

You can install `scMultiMap` from GitHub using `devtools`:

``` r
# Load devtools for installing R packages from GitHub
library(devtools)

# Install scMultiMap from GitHub
install_github("ChangSuBiostats/scMultiMap")
```

## Vignettes

The following vignettes provide detailed use cases for `scMultiMap`:

1. [Introduction to scMultiMap](https://changsubiostats.github.io/scMultiMap/articles/scMultiMap.html): 
  Learn how to infer peak-gene associations in cell types using 10x Multiome data on PBMC.

2. [scMultiMap for disease-control studies](https://changsubiostats.github.io/scMultiMap/articles/disease_control.html):
  Identify differentially associated peak-gene pairs in disease-control studies.

3. [scMultiMap for integrative analysis with GWAS results](https://changsubiostats.github.io/scMultiMap/articles/GWAS.html): 
  Integrate `scMultiMap` results with genome-wide association studies (GWAS) to explore the regulatory roles of GWAS variants in disease-associated cell types.

## scMultiMap_analysis for reproducibility

To reproduce the analysis in scMultiMap article, please visit our dedicated GitHub repository containing the source code used in the paper: [scMultiMap_analysis](https://github.com/ChangSuBiostats/scMultiMap_analysis).


## Contact us

For issues or feature requests, please visit [GitHub Issues](https://github.com/ChangSuBiostats/scMultiMap/issues). If an issue remains unanswered for a while, you are welcome to email the maintainer at chang.su@emory.edu.


## Updates

[04/26/2026] Published at Nature Communications.

[![GitHub Repo](https://img.shields.io/badge/GitHub-Repo-blue.svg)](https://github.com/ChangSuBiostats/scMultiMap)
[![DOI](https://zenodo.org/badge/926214227.svg)](https://doi.org/10.5281/zenodo.14948456)
