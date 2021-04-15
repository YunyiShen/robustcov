# RobustOmega
Robust precision matrix estimation, an R package generated from final project of STAT 760 at UW Madison in Spring 2021. Based on the review of P.-L. Loh and X. L. Tan. (2018)

There are two branches, namely `master` and `simulation`. The final project (and report) is in `simulation` branch while `master` branch is a clean R package. 

To install:

```r
devtools::install_github("YunyiShen/RobustOmega")
```

There are in total 4 robust covariance and 3 correlation estimation implemented, namely:

- `corSpearmanmat`: Spearman correlation
- `corKendallmat`: Kendall's tau
- `corQuadrantmat`: Quadrant correlation coefficients
- `covGKmat`: Gnanadesikan-Kettenring estimator by Tarr et al. (2015) and Oellerer and Croux (2015)
- `covSpearmanUmat`: SpearmanU covariance estimator by P.-L. Loh and X. L. Tan. (2018), The pairwise covariance matrix estimator proposed in Oellerer
and Croux (2015), where the MAD estimator is combined with Spearman’s
rho
- `covOGKmat`: Orthogonalized Gnanadesikan-Kettenring (OGK) estimator by Maronna, R. A. and Zamar, R. H. (2002)
- `covNPDmat`: Nearest Positive (semi)-Definite projection of the pairwise covariance matrix estimator considered in Tarr et al. (2015). 

P.-L. Loh and X. L. Tan. (2018) then used these robust estimates in Graphical Lasso (package `glasso`) or Quadratic Approximation (package `QUIC`) to obtain sparse solutions to precision matrix.


With `QUIC`, a function `robQUIC` stand for robust QUIC is implmented. It has build in cross validation described in P.-L. Loh and X. L. Tan. (2018), for instance, to use the method with cross validation:

```r
robQUIC(data=matrix(rnorm(100),20,5), covest = cov,CV=TRUE)
```

Where `data` should be a matrix and `covest` should be a function that estimate the covariance e.g. anyone mentioned above. The result list contains everything from `QUIC` output with the optimal tuning parameter found by cross validation. One can also decide fold by setting `fold` in `robQUIC`. For more details see `?robQUIC`.  

However, we have `QUIC` lost its support on CRAN, there is a branch on GitHub namely `cran/QUIC` still available to be installed. I made the package depends on the remote `QUIC` on github, so if one install the package by `devtools::install_github` then `QUIC` should not a problem.

## Simulations

I will install the package to use all the routines implemented. 

The simulations are on reconstruction of precision matrix using different robust covariance estimators fed into `glasso` and `Tiger`. 

The simulation codes are in `./simulation` while all other codes are similar to an R package. 

Several main things are considered:
**Graph Structure**

- Banded: Basically an AR(1) model
- Sparse: Sparse, with all entries positive
- Dense: Well, a dense matrix

They are all implemented in `./R/graph_generator.R`

**Generating distribution**

- Contaminated normal: normal with a normal contamination, in `./R/simu_data`
- Multivariate t: just multi-variate t, implemented in C++, `./src/t_distns.cpp`
- Alternative t: similar to multivariate t, but the scaling chi_sq distribution is different to each dimension, also in C++

**Algorithm**

As described above

## Message to repeater

I suggest not dealling with the implementations in C++ if you are not comfortable with it, several robust estimators are also implemented in base R or package `ccaPP`. I will try to claim my implementation is in C++ rather than R so you might be able to still use R. If you are a matlab user, read the documentation of C++ library `armadillo`, majority of the functions I used can be directly translated into matlab. Sorry I am not familiar to other options. 