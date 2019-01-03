
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- To add a badge  -->
<!-- [![Travis-CI Build Status](https://travis-ci.org/geanders/countyweather.svg?branch=master)](https://travis-ci.org/geanders/countyweather) -->
<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/countyweather)](https://cran.r-project.org/package=countyweather) -->
MclustSepCov
============

This is a README file of the R package *MclustSepCov*. In our paper, not yet available, we analyze the multivariate longitudinal data by applying model-based clustering. We assume the Gaussian mixture model with separable covariance structure where one of covariance factors is assumed to have a temporal structure.

Installation of the package
---------------------------

To install our package, you may simply execute the following codes:

``` r
# install.packages("devtools")
devtools::install_github("inmybrain/MclustSepCov", subdir = "MclustSepCov") # don't forget specify subdir!
```

If you come across a problem like [this](https://github.com/r-lib/remotes/issues/130), please track the issue to handle it.

Or you can install the source file using the command line after downloading it from [here](https://drive.google.com/uc?export=download&id=1l2q381dgUCr5uBN2Pp6ZjsDOYfF_U3mP);

``` bash
R CMD INSTALL MclustSepCov_1.0.tar.gz
```

A basic example of using the package
------------------------------------

We give a toy example to apply the main function `Mclust_SEP_cpp`. First, generate 40 samples from the Gaussian mixture model with two components. One of the component has zero mean vector, while the other a vector parallel to $(1, \\cdots, 1)^{\\rm T}$ with modulus 5. We assume the heterogeneous covariance structure between them.

``` r
# install.packages("mvtnorm")
library("MclustSepCov")
K <- 2
p <- 2
q <- 3
U <- lapply(1:K, function(noarg) getCovariance(p, 0.3, "AR"))
V <- lapply(1:K, function(noarg) getCovariance(q, 0.2, "CS"))
Sigma <- Map(kronecker, U, V) # separable covariance matrix
mu <- list(rep(0, p * q), 5 / sqrt(p*q) * rep(1, p * q)) # distinct mean vectors
Y <- vector(mode = "list", length = K)
set.seed(6) # to make results reproducible
for(i in 1:K){
 Y[[i]] <- mvtnorm::rmvnorm(n = 20, mean = mu[[i]], sigma = Sigma[[i]])
}
```

Each list has 20 samples from each Gaussian component. Candidate models we consider here are the Gaussian mixture distribution with the number of mixture components within {1, 2, 3} and with any of available covariance models.

``` r
fit <- Mclust_SEP_cpp(Y = Reduce(rbind, Y), p = p, q = q, Ks = 1:3, type_cov = "all", save_fit = TRUE)
```

It contains a table showing BIC values for all fitted models.

``` r
fit$bic_table
```

The best model with 2 components and covariance model 'EEE-ECS' is saved in `best_model`

``` r
fit$best_model
```

If `save_fit=TRUE` (default is `FALSE`), then one can access to all fitted models; for example,

``` r
## VVV-VUN
fit$res_Mclust_SEP$`VVV-VUN`$`1`
fit$res_Mclust_SEP$`VVV-VUN`$`2`
fit$res_Mclust_SEP$`VVV-VUN`$`3`

## EEE-EAR
fit$res_Mclust_SEP$`EEE-EAR`$`1`
fit$res_Mclust_SEP$`EEE-EAR`$`2`
fit$res_Mclust_SEP$`EEE-EAR`$`3`
```

Notes
-----

-   For available covariance structures, see the help page;

``` r
?Mclust_SEP_cpp
```

-   As for initial assignment of cluster membership, each sample is assigned randomly to clusters.

Issues
------

We are happy to troubleshoot any issue with the package;

-   please contact to the maintainer by <seongohpark6@gmail.com>, or

-   please open an issue in the github repository.

<!-- ## Error and warning messages you may get -->
<!-- ## References  -->
