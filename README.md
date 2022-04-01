
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SPCompute

<!-- badges: start -->

<!-- badges: end -->

The goal of SPCompute is to compute power and sample size for
replication GWAS study, while accommodates different kinds of covariate
effects. The methodology used in the software is described in [this
paper](https://arxiv.org/abs/2203.15641) by Ziang Zhang and Lei Sun. The
detailed implementation guideline can be found in the vignette of this
package.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AgueroZZ/SPCompute")
```

## Example

This is a basic example which shows you how to solve a common problem of
computing power for genetic association testing with a binary trait:

``` r
library(SPCompute)
## basic example code
parameters <- list(preva = 0.2, pG = 0.3, pE = 0.3, gammaG = 0.1, betaG = 0.1, betaE = 0.3)
Compute_Power(parameters, n = 8000, response = "binary", covariate = "none")
#> [1] 0.6404552
```
