// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// nearPPSD
arma::mat nearPPSD(arma::mat X, const float eigenTol, const float convTol, const float psdTol, const int maxit);
RcppExport SEXP _RobustOmega_nearPPSD(SEXP XSEXP, SEXP eigenTolSEXP, SEXP convTolSEXP, SEXP psdTolSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const float >::type eigenTol(eigenTolSEXP);
    Rcpp::traits::input_parameter< const float >::type convTol(convTolSEXP);
    Rcpp::traits::input_parameter< const float >::type psdTol(psdTolSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(nearPPSD(X, eigenTol, convTol, psdTol, maxit));
    return rcpp_result_gen;
END_RCPP
}
// covGKmat
arma::mat covGKmat(const arma::mat& data);
RcppExport SEXP _RobustOmega_covGKmat(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(covGKmat(data));
    return rcpp_result_gen;
END_RCPP
}
// corSpearmanmat
arma::mat corSpearmanmat(const arma::mat& data);
RcppExport SEXP _RobustOmega_corSpearmanmat(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(corSpearmanmat(data));
    return rcpp_result_gen;
END_RCPP
}
// corKendallmat
arma::mat corKendallmat(const arma::mat& data);
RcppExport SEXP _RobustOmega_corKendallmat(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(corKendallmat(data));
    return rcpp_result_gen;
END_RCPP
}
// corQuadrantmat
arma::mat corQuadrantmat(const arma::mat& data);
RcppExport SEXP _RobustOmega_corQuadrantmat(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(corQuadrantmat(data));
    return rcpp_result_gen;
END_RCPP
}
// covSpearmanUmat
arma::mat covSpearmanUmat(const arma::mat& data);
RcppExport SEXP _RobustOmega_covSpearmanUmat(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(covSpearmanUmat(data));
    return rcpp_result_gen;
END_RCPP
}
// covOGKmat
arma::mat covOGKmat(const arma::mat& data);
RcppExport SEXP _RobustOmega_covOGKmat(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(covOGKmat(data));
    return rcpp_result_gen;
END_RCPP
}
// covNPDmat
arma::mat covNPDmat(const arma::mat& data, const float eigenTol, const float convTol, const float psdTol, const int maxit);
RcppExport SEXP _RobustOmega_covNPDmat(SEXP dataSEXP, SEXP eigenTolSEXP, SEXP convTolSEXP, SEXP psdTolSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const float >::type eigenTol(eigenTolSEXP);
    Rcpp::traits::input_parameter< const float >::type convTol(convTolSEXP);
    Rcpp::traits::input_parameter< const float >::type psdTol(psdTolSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(covNPDmat(data, eigenTol, convTol, psdTol, maxit));
    return rcpp_result_gen;
END_RCPP
}
// raltert
arma::mat raltert(int n, const arma::mat& Omega, int nu);
RcppExport SEXP _RobustOmega_raltert(SEXP nSEXP, SEXP OmegaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(raltert(n, Omega, nu));
    return rcpp_result_gen;
END_RCPP
}
// rmvt
arma::mat rmvt(int n, const arma::mat& Omega, int nu);
RcppExport SEXP _RobustOmega_rmvt(SEXP nSEXP, SEXP OmegaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvt(n, Omega, nu));
    return rcpp_result_gen;
END_RCPP
}
// rmvnorm
arma::mat rmvnorm(int n, const arma::mat& Omega);
RcppExport SEXP _RobustOmega_rmvnorm(SEXP nSEXP, SEXP OmegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvnorm(n, Omega));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RobustOmega_nearPPSD", (DL_FUNC) &_RobustOmega_nearPPSD, 5},
    {"_RobustOmega_covGKmat", (DL_FUNC) &_RobustOmega_covGKmat, 1},
    {"_RobustOmega_corSpearmanmat", (DL_FUNC) &_RobustOmega_corSpearmanmat, 1},
    {"_RobustOmega_corKendallmat", (DL_FUNC) &_RobustOmega_corKendallmat, 1},
    {"_RobustOmega_corQuadrantmat", (DL_FUNC) &_RobustOmega_corQuadrantmat, 1},
    {"_RobustOmega_covSpearmanUmat", (DL_FUNC) &_RobustOmega_covSpearmanUmat, 1},
    {"_RobustOmega_covOGKmat", (DL_FUNC) &_RobustOmega_covOGKmat, 1},
    {"_RobustOmega_covNPDmat", (DL_FUNC) &_RobustOmega_covNPDmat, 5},
    {"_RobustOmega_raltert", (DL_FUNC) &_RobustOmega_raltert, 3},
    {"_RobustOmega_rmvt", (DL_FUNC) &_RobustOmega_rmvt, 3},
    {"_RobustOmega_rmvnorm", (DL_FUNC) &_RobustOmega_rmvnorm, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_RobustOmega(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}