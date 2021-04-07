#ifndef NPD_PROJ_H
#define NPD_PROJ_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
mat ProjPattern(mat X, const mat X0, const vec un);
mat ProjPSD(const mat R, const int n, const float eigenTol);
mat nearPPSD(mat X, const float eigenTol = 1e-06, const float convTol = 1e-07, 
             const float psdTol = 1e-08, const int maxit = 1000);
#endif