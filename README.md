[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# RICE
This repository contains an R package for creating polygenic risk scores from rare variants using RICE. Specifically, this package creates a polygenic risk score using RICE-RV (as denoted in the manuscript) with significant rare variant sets and their burden scores extracted from an AGDS file. A detailed example of this is provided in the <a href="https://jwilliams10.github.io/RICE_Vignette">vignette</a>. **Note**, if an individual has burden scores from a source other than an AGDS file, a rare variant PRS following the RICE-RV steps can still be created by starting at the Penalized Regression step in the vignette.

## RICE Overview
![alt text](https://github.com/jwilliams10/RICE/blob/main/docs/RICE_Structure.jpg?raw=true)

## Prerequisites
**RICE** can be viewed as an extension to <a href="https://github.com/xihaoli/STAARpipeline">**STAARpipeline**</a>. We suggest users familiarized with STAARpipeline and it's <a href="https://github.com/xihaoli/STAARpipeline-Tutorial">tutorial</a> before constructing rare variant PRSs using RICE.

<a href="https://www.r-project.org">R</a> (recommended version >= 4.0.0)

**RICE** imports R packages Rcpp, SCANG, dplyr, SeqArray, SeqVarTools, GenomicFeatures, TxDb.Hsapiens.UCSC.hg38.knownGene, GMMAT, GENESIS, Matrix, methods, caret, glmnet, SuperLearner, STAAR, and STAARpipeline. These dependencies should be installed before installing **RICE**.

## Installation
```R
library(devtools)
devtools::install_github("jwilliams10/RICE")
```

Installs in approximately 10 seconds on a Macbook Pro M1 if prerequisites are already installed.

## Usage
Please see the <a href="docs/RICE_0.1.0.pdf">**RICE** user manual</a> for detailed usage of RICE package.

An uncompiled <a href="vignettes/RICE_Vignette.Rmd">vignette</a> and a compiled <a href="https://jwilliams10.github.io/RICE_Vignette">vignette</a> is available that showcases how to create a rare variant PRS with part of chromosome 22 for 1000 Genome data. Data for the example in the vignette is downloadable through <a href="https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/Y7AXYA">Harvard Dataverse</a>. The vignette takes approximately 70 seconds to compile on a Macbook Pro M1.

## Support
Please direct any problems or questions to Jacob Williams <jacob.williams@nih.gov>.

## Citation
If you use **RICE** in your work, please cite:

Jacob Williams*, Tony Chen, Xing Hua, Wendy Wong, Kai Yu, Peter Kraft, Xihao Li*, & Haoyu Zhang*. (2024). **Integrating Common and Rare Variants Improves Polygenic Risk Prediction Across Diverse Populations**. _medRxiv_. DOI: <a href="https://doi.org/10.1101/2024.11.05.24316779">2024.11.05.24316779</a> 

## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
