#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


// ' Alternative multivariate t-distribution
// ' This routine samples alternative multivarate t-distribution
// ' @param n sample size
// ' @param Omega precision matrix of dimension p by p (not covariance)
// ' @param nu degree of freedom
// ' @return a matrix with dimension n by p, each row is a sample
// ' @export
// [[Rcpp::export]]
arma::mat raltert(int n, const arma::mat & Omega, int nu){
	mat L = chol(Omega); 
    int p = Omega.n_cols;
    mat res(p,n, arma::fill::randn);
	arma::solve(res, L, res); // normal (0, sigma)
	mat temp_chi = arma::chi2rnd( nu, arma::size(res) );
    res = res / (arma::sqrt(temp_chi/nu));
    return(res.t());
}


// ' Multivariate t-distribution
// ' This routine samples multivarate t-distribution
// ' @param n sample size
// ' @param Omega precision matrix of dimension p by p (not covariance)
// ' @param nu degree of freedom
// ' @return a matrix with dimension n by p, each row is a sample
// ' @export
// [[Rcpp::export]]
arma::mat rmvt(int n, const arma::mat & Omega, int nu){
	mat L = chol(Omega); 
    int p = Omega.n_cols;
    mat res(p, n, arma::fill::randn);
	arma::solve(res, L, res); // normal (0, sigma)
	vec temp_chi = chi2rnd( nu, n );
    res.each_row() /= arma::sqrt(temp_chi.t()/nu);
    return(res.t());
}


// ' Multivariate normal distribution with 0 mean
// ' This routine samples multivarate normal distribution of mean 0 from precision matrix
// ' @param n sample size
// ' @param Omega precision matrix of dimension p by p (not covariance)
// ' @return a matrix with dimension n by p, each row is a sample
// ' @export
// [[Rcpp::export]]
arma::mat rmvnorm(int n, const arma::mat & Omega){
	mat L = chol(Omega); 
    int p = Omega.n_cols;
    mat res(p, n, arma::fill::randn);
    arma::solve(res, L, res);
    return(res.t());
}