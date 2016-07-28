#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

/*Compute loglikelihood of ghk */
// [[Rcpp::export]] 
List ghksimM(arma::vec mu, arma::mat R, arma::vec lower, arma::vec upper, int M) {
  int dim = R.n_cols;double z_[dim];
  arma::mat L = arma::chol(R, "lower");
  arma::mat sim(dim, M); 
  arma::vec logweight(M);
  double  MU_,  eta_, gamma_, ans_;
  const double EPS = std::numeric_limits<double>::min(); 
  double* eta = & eta_;  double* ans = &ans_; 
  double* gamma = & gamma_; double* MU = & MU_; double* z = &z_[0];
  arma::vec weight(dim);
  for(int i=0; i< M; i++){
     for(int j=0;j<dim;j++) {       
    *MU = 0; 
    for(int k=0;k<j;k++)  *MU += L(j,k)*(*(z+k));
    *eta = Rcpp::stats::pnorm_1((lower(j)-*MU - mu(j))/L(j,j),0,1,0);
    *gamma = Rcpp::stats::pnorm_1((upper(j)-*MU - mu(j))/L(j,j),0,1,0);
    double u = Rf_runif(0, 1);
    *ans=u*(*gamma)+(1-u)*(*eta);
    if(*ans > 1-EPS) *ans= 1-EPS;
    if(*ans < EPS) *ans= EPS;
    *(z+j) = Rcpp::stats::qnorm_1(*ans,0,1,0);
    sim(j, i) = mu(j)+ *MU + L(j,j)*(*(z+j));
    weight(j) = log(*gamma - *eta);
  }
  logweight(i) =  accu(weight);
  }
  return List::create(Named("simulation") = sim, Named("logweight") = logweight);
}


