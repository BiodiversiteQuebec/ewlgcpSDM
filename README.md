# ewlgcpSDM

ewlgcpSDM (**e**ffort-**w**eighted **L**og-**G**aussian **C**ox **P**rocess **s**pecies **d**istribution **m**odels) is a package for inferring species distributions using presence-only data and point process models. In particular, it implements the Log-Gaussian Cox process and adjusts for sampling effort by weighting the thinned surface through and estimate of the effort. The effort surface is derived through a method analogous to Target-Group Background selection [Phillips *et al*. 2009](https://doi.org/10.1890/07-2153.1), whereas the sampling effort is approximated by grouping all observations of a target group of similar species. The implementation is adapted from the method proposed by [Simpson *et al*. (2016)](https://doi.org/10.1093/biomet/asv064) using [INLA](https://www.r-inla.org/) ([Rue *et al*. 2009](https://doi.org/10.1111/j.1467-9868.2008.00700.x)).


## Installation

First, INLA (not on CRAN) has to be installed (see [here](https://www.r-inla.org/download-install) for instructions) or do:

```r
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
```

Then, ewlgcpSDM can be installed directly from GitHub like this:
```r
remotes::install_github("frousseu/ewlgcpSDM")
```

## References

Phillips, S. J., Dudík, M., Elith, J., Graham, C. H., Lehmann, A., Leathwick, J. and Ferrier, S. 2009. Sample selection bias and presence-only distribution models: implications for background and pseudo-absence data. Ecological Applications, 19(1): 181-197 [https://doi.org/10.1890/07-2153.1](https://doi.org/10.1890/07-2153.1)

Rue, H., Martino, S. and Chopin, N. 2009. Approximate Bayesian Inference for Latent Gaussian models by using Integrated Nested Laplace Approximations, Journal of the Royal Statistical Society Series B: Statistical Methodology, 71(2): 319–392 [https://doi.org/10.1111/j.1467-9868.2008.00700.x](https://doi.org/10.1111/j.1467-9868.2008.00700.x)

Simpson, D., Illian, J. B., Lindgren, F., Sørbye, S. H. and Rue, H. 2016. Going off grid: computationally efficient inference for log-Gaussian Cox processes. Biometrika, 103(1): 49-70 [https://doi.org/10.1093/biomet/asv064](https://doi.org/10.1093/biomet/asv064)
