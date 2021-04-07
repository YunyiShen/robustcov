#ifndef FASTCORKENDALL_H
#define FASTCORKENDALL_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
//#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// functions to be used within C++
double fastCorKendall(const vec& x, const vec& y, const uword & n);

#endif
