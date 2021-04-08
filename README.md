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

P.-L. Loh and X. L. Tan. (2018) then used these robust estimates in Graphical Lasso (package `glasso`) or Quadratic Approximation (package `QUIC`) to obtain sparse solutions to precision matrix