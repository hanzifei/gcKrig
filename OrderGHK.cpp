#include <Rmath.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace arma; 
using namespace Rcpp;


// [[Rcpp::export]]
double ghkLiko(arma::vec mu, arma::mat R, arma::vec a, arma::vec b, int nrep){
  
 const double eps = std::numeric_limits<double>::epsilon();  
 const double EPS = std::numeric_limits<double>::min(); 
 int n = R.n_cols;
 arma::mat L = R; arma::vec ap = a; arma::vec bp = b; 
 arma::vec y;  y.zeros(n);
 double cii, ckk, dem, am, bm, de, ai, bi,s;
 int im; 
 arma::mat tmp; arma::mat mattmp;

for(int k=0; k<n; k++){
  im = k;  ckk = 0;  dem = 1; s = 0;
  
  for(int i=k; i<n; i++){
    
    if(L(i,i) > eps){
      cii = sqrt(R::fmax2(L(i,i), 0 ));
      if(k == 0){s = 0;}
      if(i > 0 && k > 0){ 
        mattmp = L(i, span(0,k-1))*y(span(0,k-1)); s = mattmp(0,0);
      }
      
      ai = (ap(i)-s)/cii;  bi = (bp(i)-s)/cii; 
      de = R::pnorm(bi,0,1,1,0) - R::pnorm(ai,0,1,1,0);
      
      if(de<=dem){ ckk = cii;  dem = de;  am = ai;  bm = bi; im = i;}
    }
  }
  
  if(im > k){
    ap.swap_rows(im,k); bp.swap_rows(im,k);
    L(im, im) = L(k,k);
    
    if(k > 0){
      L(span::all,span(0,k-1)).swap_rows(im,k);
    }
    
    if( n-im >= 2){
      L(span(im+1,n-1),span::all).swap_cols(im,k);
    }
    
    if(im - k >= 2){ 
      tmp = L(span(k+1,im-1),k); L(span(k+1,im-1),k) = L(im,span(k+1,im-1)).t();
      L(im,span(k+1,im-1)) = tmp.t();
    }
  }
  
  if(n - k >= 2)  L(k,span(k+1,n-1)).zeros();
  
  if(ckk > eps*k){
    L(k,k) = ckk; 
    for(int i=k+1; i<n; i++){
      L(i,k) = L(i,k)/ckk;
      L(i, span(k+1,i)) = L(i, span(k+1,i))-L(i,k)*L(span(k+1,i),k).t();
    }
    if(std::abs(dem) > eps){
      y(k) = (R::dnorm(am,0,1,0) - R::dnorm(bm,0,1,0) )/ (dem);
    }else{
      if(am < -10) {y(k) = bm;}
      else if(bm > 10) {y(k) = am;}
      else{y(k) = (am+bm)/2;}
    }
  }
  else{L(span(k,n-1), k).zeros(); y(k) = 0; }
  
  }

//  return Rcpp::List::create(Rcpp::Named("L") = L, 
//                            Rcpp::Named("a") = ap,
//                            Rcpp::Named("b") = bp);  
    
  double z_[n];
  double res_ = 0, prod, MU_, eta_, gamma_, ans_;

  double* eta = & eta_; double* res = & res_; double* ans = &ans_;
  double* gamma = & gamma_; double* MU = & MU_; double* z = &z_[0];
  
  for(int ii=0; ii<nrep; ii++) {
    
    prod = 1.0;
    
    for(int jj=0; jj<n; jj++) {    
      
      *MU=0; 
      
      for(int kk=0; kk<jj; kk++)  *MU += L(jj,kk)*(*(z+kk));
      
      *eta = Rcpp::stats::pnorm_1((ap(jj)-*MU-mu(jj))/L(jj,jj),0,1,0);
      
      *gamma = Rcpp::stats::pnorm_1((bp(jj)-*MU-mu(jj))/L(jj,jj),0,1,0);
      
      prod *= *gamma-*eta; 
      double u = Rf_runif(0, 1);
      *ans=u*(*gamma)+(1-u)*(*eta);
      if(*ans > 1-EPS) *ans= 1-EPS;
      if(*ans < EPS) *ans= EPS ;
      *(z+jj) = Rcpp::stats::qnorm_1(*ans,0,1,0);
    }
    *res += prod;
  }
  *res /= nrep;
  if(*res == 0){
    return -1e6;
  }
  else{
    return log(*res);
  }
}
  
  
  
  
  

  
  
  
  
  
  
  
  
