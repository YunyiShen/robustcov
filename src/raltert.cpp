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
    mat res(n,p); 
    for(int i = 0 ; i < n ; i++){
        vec temp_n(p, arma::fill::randn);
	    temp_n = arma::solve(L, temp_n); // normal (0, sigma)
	    vec temp_chi = chi2rnd( nu, p );
        res.row(i) = trans( temp_n / (arma::sqrt(temp_chi/nu)));

    }
    return(res);
}