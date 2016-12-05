# BridgeR2
[![Travis-CI Build Status](https://travis-ci.org/Imamachi-n/BridgeR2.svg?branch=master)](https://travis-ci.org/Imamachi-n/BridgeR2)
[![codecov](https://codecov.io/gh/Imamachi-n/BridgeR2/branch/master/graph/badge.svg)](https://codecov.io/gh/Imamachi-n/BridgeR2)  
  
## Overview
BRIC-seq is a genome-wide approach for determining RNA stalibity in mammalian cells. `bridger2` provides a series of functions for performing a comprehensive BRIC-seq data analysis. After estimating the RPKM values for all genes from your BRIC-seq fastq files, you can easily analyze your BRIC-seq data using bridger2 R package. 

To make that happen, `bridger2`:
* Checks the quality of your BRIC-seq data.

* Normalizes RPKM values of your BRIC-seq data.

* Calculates RNA half-life for each transcript

* Compares RNA half-lives between two conditions.

* Displays RNA decay curve using a web browser (powered by shiny).

## Installation
```r
# The the development version from GitHub:
# install.packages("devtools")
devtools::install_github("Imamachi-n/BridgeR2")
```

## Quick start
Here I show the most basic step for analyzing your BRIC-seq data. This step require `matrix` object (named `RNA_halflife_comparison` in this case) of the RPKM values from your BRIC-seq data. `BridgeRCore` function returns `data.table` object including RNA half-life, R2 and the selected fitting model.
```r
halflife_table <- BridgeRCore(RNA_halflife_comparison)
```
