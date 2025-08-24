# deepMatchR

Tools for HLA Testing and Matching

<!-- badges: start -->
[![R-CMD-check](https://github.com/BorchLab/deepMatchR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BorchLab/deepMatchR/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/BorchLab/deepMatchR/graph/badge.svg)](https://app.codecov.io/gh/BorchLab/deepMatchR)
<!-- badges: end -->

## Introduction

<img align="right" src="https://github.com/BorchLab/deepMatchR/blob/main/www/deepMatchR_hex.png" width="305" height="352">

There are currently several computational approaches to quantifying the risk of the development of donor-specific antibodies during organ transplantation. These include [HLAmatchmaker](http.www.epitopes.net/) for eplet quantification and [PIRCH-II](https://www.pirche.com/) for CD4+ T cell epitope prediction, which have demonstrated predictive ability across the literature. Newer deep learning methods for structure predictions, eplet/epitope immunogenicity estimates, and classification can be leveraged to produce a clinical tool for patients. deepMatchR aims to be a centralized repository for tools and models to help in assisting HLA matching. It now includes the `spiDeconvolute` function for reconciling single antigen bead (SAB) and panel reactive antibody (PRA) results.

## System requirements 

deepMatchR has been tested on R versions >= 4.0. Please consult the DESCRIPTION file for more details on required R packages. deepMatchR has been tested on OS X and Windows platforms.

## Installation

To run deepMatchR, open R and install deepMatchR from github: 

```r
devtools::install_github("ncborcherding/deepMatchR")
```
***
## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/ncborcherding/deepMatchR/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 
Alternatively, an example with the internal **deepMatchR_example** would 
be extremely helpful.

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/ncborcherding/deepMatchR/issues).

#### [Pull Requests](https://github.com/ncborcherding/deepMatchR/pulls) are welcome for bug fixes, new features, or enhancements.

