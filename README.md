# deepMatchR
Deep learning-based approaches to HLA matching

## Introduction

<img align="right" src="https://github.com/BorchLab/deepMatchR/blob/main/www/deepMatchR_hex.png" width="305" height="352">

There are currently several computational approaches to quantifying the risk of the development of donor-specific antibodies during organ transplantation. These include [HLAmatchmaker](http.www.epitopes.net/) for eplet quantification and [PIRCH-II](https://www.pirche.com/) for CD4+ T cell epitope prediction, which have demonstrated predictive ability across the literature. Newer deep learning methods for structure predictions, eplet/epitope immunogenicity estimates, and classification can be leveraged to produce a clinical tool for patients. deepMatchR aims to be a centralized repository for tools and models to help in assisting HLA matching. It now includes the `spiDeconvolute` function for reconciling single antigen bead (SAB) and panel reactive antibody (PRA) results.

## System requirements 

deepMatchR has been tested on R versions >= 4.0. Please consult the DESCRIPTION file for more details on required R packages. deepMatchR has been tested on OS X and Windows platforms.

**keras** is necessary to use the autoencoder function (this includes the set up of the tensorflow environment in R):

```r
##Install keras
install.packages("keras")

##Setting up Tensor Flow
library(reticulate)
conda_create("r-reticulate") ##If first time using reticulate
use_condaenv(condaenv = "r-reticulate", required = TRUE)
library(tensorflow)
install_tensorflow()
```

An alternative to this approach above (especially if you want to avoid conda) is to use reticulate to generate a virtualenv, using ```virtualenv_create()``` and subsequently installing the above python packages using ```virtualenv_install()```.

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

