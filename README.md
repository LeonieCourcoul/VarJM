
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VarJM

<!-- badges: start -->
<!-- badges: end -->

The goal of VarJM is to estimate joint model with subject-specific
variability.

The global function is FlexVar_JM. It handles to estimate joint model
with a marker which has a subject-specific variability and competing
events with the possibility to take into account the left truncation.

## Installation

You can install the development version of VarJM from
[GitHub](https://github.com/) with:

``` r
#devtools::install_github("LeonieCourcoul/VarJM")
```

## Exemple

This exemple is performed with a simulated toy dataset with two
competing risks and a Weibull baseline risk function.

``` r
#exemple <- FlexVar_JM(formFixed = y~visit,
#                      formRandom = ~ visit,
#                      formGroup = ~ID,
#                      formSurv = Surv(time, event ==1 ) ~ 1,
#                      timeVar = "visit",
#                      nb.e.a = 2,
#                      data.long = Data_exemple,
#                      variability_hetero = TRUE,
#                      sharedtype = "CV",
#                      hazard_baseline = "Weibull",
#                      competing_risk = TRUE,
#                      formSurv_CR = Surv(time, event ==2 ) ~ 1,
#                      hazard_baseline_CR = "Weibull",
#                      sharedtype_CR = "CV",
#                      S1 = 1000,
#                      S2 = 8000,
#                      nproc = 5
#                      )
#
#exemple$table.res
```
