[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# RICE
This an R package for creating polygenic risk scores from rare variants using RICE. Creates a rare variant polygenic risk score using burden scores extracted from a ADGS file. Burden scores are extracted for specific rare variant sets using the Burden_Score function. A rare variant PRS is constructed using the RareVariantPRS function.

## RICE Overview
![alt text](https://github.com/jwilliams10/RICE/blob/main/docs/RICE_Structure.jpg?raw=true)

## Prerequisites
**RICE** can be viewed as an extension to **STAARpipeline** (https://github.com/xihaoli/STAARpipeline). We suggest users familiarized with STAARpipeline and it's tutorial (https://github.com/xihaoli/STAARpipeline-Tutorial) before constructing rare variant PRSs using RICE.

<a href="https://www.r-project.org">R</a> (recommended version >= 4.0.0)

**RICE** imports R packages Rcpp, SCANG, dplyr, SeqArray, SeqVarTools, GenomicFeatures, TxDb.Hsapiens.UCSC.hg38.knownGene, GMMAT, GENESIS, Matrix, methods, caret, glmnet, SuperLearner, STAAR, and STAARpipeline. These dependencies should be installed before installing **RICE**.

## Installation
```R
library(devtools)
devtools::install_github("jwilliams10/RICE")
```

## Usage

## Support
Please direct any problems or questions to Jacob Williams <jacob.williams@nih.gov>.

## Citation

If you use **RICE** in your work, please cite:

...

## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
