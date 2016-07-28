#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace arma;

//Come from RcppArmadillo example for computing log-density
// [[Rcpp::export]]

double ldgc(arma::rowvec u, arma::mat sigma) { 
  
  int dimn = sigma.n_rows; 
  
  arma::mat RI = arma::inv(sigma)-arma::speye(dimn, dimn);
  
  arma::rowvec phi(dimn);
  
  for(int i=0; i<dimn; i++){ 
    phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);
    }
  
  arma::mat outtmp = phi*RI*phi.t();
    
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
    
  double out = -0.5*(logdet+outtmp(0));
   
  return(out);
}
 