#include <Rmath.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace arma; 
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List orderGHK(arma::mat R, arma::vec a, arma::vec b){
  
  double eps = std::numeric_limits<double>::epsilon();   int n = R.n_cols;
  arma::mat c = R; arma::vec ap = a; arma::vec bp = b; 
  arma::vec y;  y.zeros(n);
  double cii, ckk, dem, am, bm, de, ai, bi,s;
  int im;
  arma::mat tmp;
  arma::mat mattmp;
  

  for(int k=0; k<n; k++){
      im = k;  ckk = 0;  dem = 1; s = 0;
    
    for(int i=k; i<n; i++){
      
      if(c(i,i) > eps){
        cii = sqrt(R::fmax2(c(i,i), 0 ));
        
      if(k == 0){s = 0;}
      if(i > 0 && k > 0){ 
          mattmp = c(i, span(0,k-1))*y(span(0,k-1)); s = mattmp(0,0);
            }
        
    ai = (ap(i)-s)/cii;  bi = (bp(i)-s)/cii; 
    de = R::pnorm(bi,0,1,1,0) - R::pnorm(ai,0,1,1,0);
       
    if(de<=dem){ ckk = cii;  dem = de;  am = ai;  bm = bi; im = i;}
      }
    }
    
    if(im > k){
      ap.swap_rows(im,k); bp.swap_rows(im,k);
      c(im, im) = c(k,k);
      
      if(k > 0){
      c(span::all,span(0,k-1)).swap_rows(im,k);
      }
      
      if( n-im >= 2){
      c(span(im+1,n-1),span::all).swap_cols(im,k);
      }
     
      if(im - k >= 2){ 
      tmp = c(span(k+1,im-1),k); c(span(k+1,im-1),k) = c(im,span(k+1,im-1)).t();
      c(im,span(k+1,im-1)) = tmp.t();
     }
  }
  
  if(n - k >= 2)  c(k,span(k+1,n-1)).zeros();

    if(ckk > eps*k){
      c(k,k) = ckk; 
       for(int i=k+1; i<n; i++){
         c(i,k) = c(i,k)/ckk;
         c(i, span(k+1,i)) = c(i, span(k+1,i))-c(i,k)*c(span(k+1,i),k).t();
       }
       if(std::abs(dem) > eps){
         y(k) = (R::dnorm(am,0,1,0) - R::dnorm(bm,0,1,0) )/ (dem);
       }else{
         if(am < -10) {y(k) = bm;}
         else if(bm > 10) {y(k) = am;}
         else{y(k) = (am+bm)/2;}
       }
    }
    else{c(span(k,n-1), k).zeros(); y(k) = 0; }
    
  }
    
  return Rcpp::List::create(Rcpp::Named("L") = c, 
                            Rcpp::Named("a") = ap,
                            Rcpp::Named("b") = bp);  
}
  
  
  
  
