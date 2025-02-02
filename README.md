# scMultiMap

`scMultiMap` is an R package for inferring cell-type-specific peak-gene associations using single-cell multimodal data. It implements the statistical method proposed in th manuscript [Cell-type-specific mapping of enhancer and target genes from single-cell multimodal data](https://www.biorxiv.org/content/10.1101/2024.09.24.614814v1), currently under revision at *Nature Communications*.

To reproduce the analysis in  [Cell-type-specific mapping of enhancer and target genes from single-cell multimodal data](https://www.biorxiv.org/content/10.1101/2024.09.24.614814v1), please visit our dedicated GitHub repository containing the source code used in the paper: [scMultiMap_analysis](https://changsubiostats.github.io/scMultiMap).

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


## Contact us

For questions or feedback, please contact:
[Chang Su](www.changsu.org), <chang.su@emory.edu>

## Reference and Updates

Chang Su, Dongsoo Lee, Peng Jin and Jingfei Zhang. (2024). [Cell-type-specific mapping of enhancer and target genes from single-cell multimodal data](https://www.biorxiv.org/content/10.1101/2024.09.24.614814v1). Manuscript under revision at *Nature Communications*.
