#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace arma;

// [[Rcpp::export]]
double gqtlikNB(arma::vec count, arma::vec mu, double od, arma::mat sigma) { 
  
  int dimn = sigma.n_rows; 
  double SZ = 1/od; arma::vec PB = SZ/(SZ+mu);
  arma::mat RI = arma::inv(sigma)-arma::speye(dimn, dimn);
  arma::rowvec phi(dimn); arma::rowvec u(dimn);
  
  double loglik = 0;  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  
  for(int i=0; i<dimn; i++){ 
    u[i] = (R::pnbinom(count[i], SZ, PB[i], 1, 0) + R::pnbinom(count[i]-1, SZ, PB[i], 1, 0))/2;
    phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);
    loglik += R::dnbinom(count[i], SZ, PB[i], 1);
  }
  arma::mat outtmp = phi*RI*phi.t();
  double out = -0.5*(logdet+outtmp(0))+loglik;
  if(std::isinf(out) == true || std::isnan(out) == true){
    out = -9.532493e+14;
  }
  return(out);
}



// [[Rcpp::export]]
double gqtlikPois(arma::vec count, arma::vec mu, arma::mat sigma) { 
  int dimn = sigma.n_rows; 
  arma::mat RI = arma::inv(sigma)-arma::speye(dimn, dimn);
  arma::rowvec phi(dimn); arma::rowvec u(dimn);
  double loglik = 0;  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  for(int i=0; i<dimn; i++){ 
    u[i] = (R::ppois(count[i], mu[i], 1, 0) + R::ppois(count[i]-1, mu[i], 1, 0))/2;
    phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);
    loglik += R::dpois(count[i], mu[i], 1);
  }
  arma::mat outtmp = phi*RI*phi.t();
  double out = -0.5*(logdet+outtmp(0))+loglik;
  if(std::isinf(out) == true || std::isnan(out) == true){
    out = -9.532493e+14;
  }
  return(out);
}

double pZIPinner(double count, double mu, double od){
  double cdf = 0;
  if(count == 0){
    cdf = od/(1+od) + exp(-(1+od)*mu)/(1+od);
  }
  if(count > 0){
    cdf = R::ppois(count, (1+od)*mu, 1, 0)/(1+od) + od/(1+od);
  }
  return cdf;
}

double logdZIPinner(double count, double mu, double od){
  double pdf = 0;
  if(count == 0){
    pdf = od/(1+od) + exp(-(1+od)*mu)/(1+od);
  }
  if(count > 0){
    pdf = R::dpois(count, (1+od)*mu, 0)/(1+od);
  }
  return log(pdf);
}



// [[Rcpp::export]]
double gqtlikZIP(arma::vec count, arma::vec mu, double od, arma::mat sigma) { 
  
  int dimn = sigma.n_rows; 
  arma::mat RI = arma::inv(sigma)-arma::speye(dimn, dimn);
  arma::rowvec phi(dimn); arma::rowvec u(dimn);
  
  double loglik = 0;  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  double tmp = 0;
  
  for(int i=0; i<dimn; i++){ 
    
    if(count[i] == 0){
      u[i] = (od/(1+od) + exp(-(1+od)*mu[i])/(1+od))/2;
      tmp = log(od/(1+od) + exp(-(1+od)*mu[i])/(1+od));
    }
    if(count[i] == 1){
      u[i] = (exp(-(1+od)*mu[i])/(1+od)+R::ppois(1, (1+od)*mu[i], 1, 0)/(1+od))/2 + od/(1+od);
      tmp = R::dpois(1, (1+od)*mu[i], 1) - log(1+od);
    }
    if(count[i] > 1){
      u[i] = (R::ppois(count[i], (1+od)*mu[i], 1, 0)/(1+od)+
              R::ppois(count[i]-1, (1+od)*mu[i], 1, 0)/(1+od))/2+od/(1+od);
      tmp = R::dpois(count[i], (1+od)*mu[i], 1) - log(1+od);
    }
    phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);

    loglik += tmp;
  }
  arma::mat outtmp = phi*RI*phi.t();
  double out = -0.5*(logdet+outtmp(0))+loglik;
  if(std::isinf(out) == true || std::isnan(out) == true){
    out = -9.532493e+14;
  }
  return(out);
}



// [[Rcpp::export]]
double gqtlikBin(arma::vec count, arma::vec p, double N, arma::mat sigma) { 

  int dimn = sigma.n_rows; 
//  double SZ = 1/od; arma::vec PB = SZ/(SZ+mu);
  arma::mat RI = arma::inv(sigma)-arma::speye(dimn, dimn);
  arma::rowvec phi(dimn); arma::rowvec u(dimn);
  
  double loglik = 0;  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  
  for(int i=0; i<dimn; i++){ 
    u[i] = (R::pbinom(count[i], N, p[i], 1, 0) + R::pbinom(count[i]-1, N, p[i], 1, 0))/2;
    phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);
    loglik += R::dbinom(count[i], N, p[i], 1);
  }
  arma::mat outtmp = phi*RI*phi.t();
  double out = -0.5*(logdet+outtmp(0))+loglik;
  if(std::isinf(out) == true || std::isnan(out) == true){
    out = -9.532493e+14;
  }
  return(out);
}