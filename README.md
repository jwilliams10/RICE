# RICE

This an R package for creating polygenic risk scores from rare variants using RICE. 

## RICE Overview

![alt text](https://github.com/jwilliams10/RICE/blob/main/docs/RICE_Structure.jpg?raw=true)

## Description

Creates a rare variant polygenic risk score using burden scores extracted from a ADGS file. Burden scores are extracted for specific rare variant sets using the Burden_Score function. A rare variant PRS is constructed using the RareVariantPRS function.

## Installation

```R
library(devtools)
install_github("jwilliams10/RICE")
```
